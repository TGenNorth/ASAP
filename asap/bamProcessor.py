#!/usr/bin/env python3
# encoding: utf-8
'''
asap.bamProcessor -- Process BAM alignment files with an AssayInfo JSON file and generate XML for the results

asap.bamProcessor 

@author:     Darrin Lemmer

@copyright:  2015,2019 TGen North. All rights reserved.

@license:    ACADEMIC AND RESEARCH LICENSE -- see ../LICENSE

@contact:    dlemmer@tgen.org
'''

import sys
import os
import re
import argparse
import logging
import math

import pysam
from collections import Counter
from xml.etree import ElementTree
from skbio import DNA

from asap import dispatcher
from asap import assayInfo
from asap import __version__

__all__ = []
__updated__ = '2019-01-15'
__date__ = '2015-07-16'

DEBUG = 1
TESTRUN = 0
PROFILE = 0

low_level_cutoff = 0.01
high_level_cutoff = 0.50

def pairwise(iterable):
    from itertools import tee
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

def _write_parameters(node, data):
    for k, v in data.items():
        subnode = ElementTree.SubElement(node, k)
        subnode.text = str(v)
    return node

def _process_pileup_SMOR(pileup, amplicon, depth, proportion, mutdepth, offset, wholegenome):
    from operator import attrgetter
    pileup_dict = {}
    snp_dict = _create_snp_dict(amplicon)
    consensus_seq = ""
    snp_list = []
    breadth_positions = 0
    avg_depth_total = avg_depth_positions = 0
    amplicon_length = len(amplicon.sequence)
    depth_array = [0] * amplicon_length
    discard_array = [0] * amplicon_length
    prop_array = ["0"] * amplicon_length
    for pileupcolumn in pileup:
        base_counter = Counter()
        reads = iter(sorted(pileupcolumn.pileups, key=attrgetter('alignment.query_name')))
        #reads = iter(pileupcolumn.pileups)
        for read, pair in pairwise(reads):
            if read.alignment.query_name != pair.alignment.query_name:
                continue
            alignment = read.alignment
            if pair.is_del and read.is_del:
                base_counter.update("_") # XSLT doesn't like '-' as an attribute name, have to use '_'
            elif read.query_position and pair.query_position:
                if pair.alignment.query_sequence[pair.query_position] == alignment.query_sequence[read.query_position]:
                    depth_array[pileupcolumn.pos] += 1
                    base_counter.update(alignment.query_sequence[read.query_position])
                else:
                    discard_array[pileupcolumn.pos] += 1
                #    print("Pair does not match at reference position %d -> (%s != %s)" % (pileupcolumn.reference_pos, alignment.query_sequence[read.query_position], pair.alignment.query_sequence[pair.query_position]))
        column_depth = depth_array[pileupcolumn.pos]
        if column_depth == 0:
            consensus_seq += "_" #TODO Temporary, remove this line!
            continue
        position = pileupcolumn.pos+1
        depth_passed = False
        if column_depth > 0: #TODO: This is going to end up being specific to these TB assays, maybe have a clever way to make this line optional
            avg_depth_positions += 1
            avg_depth_total += column_depth
        if column_depth >= depth:
            breadth_positions += 1
            depth_passed = True
        ordered_list = base_counter.most_common()
        alignment_call = ordered_list[0][0]
        alignment_call_proportion = ordered_list[0][1] / column_depth
        prop_array[pileupcolumn.pos] = "%.3f" % alignment_call_proportion
        reference_call = amplicon.sequence[pileupcolumn.pos]
        if reference_call == '-':
            reference_call = '_' #Need to use '_' instead of '-' for gaps because of XSLT
        if alignment_call != reference_call:
            snp_call = alignment_call
            snp_count = ordered_list[0][1]
            snp_call_proportion = alignment_call_proportion
        elif len(ordered_list) > 1:
            snp_call = ordered_list[1][0]
            snp_count = ordered_list[1][1]
            snp_call_proportion = ordered_list[1][1] / column_depth
        else:
            snp_call = snp_count = snp_call_proportion = None
        consensus_seq += alignment_call if alignment_call_proportion >= proportion else "N"
        (proportion, low_level_cutoff, high_level_cutoff) = _compute_thresholds_SMOR(column_depth)
        translated = offset+position
        if position in snp_dict:
            for (name, reference, variant, significance) in snp_dict[position]:
                snp = {'name':name, 'position':str(translated), 'depth':str(column_depth), 'reference':reference, 'variant':variant, 'basecalls':base_counter}
                variant_proportion = base_counter[variant]/column_depth
                variant_count = base_counter[variant]
                if variant_proportion >= proportion and variant_count >= mutdepth:
                    snp['significance'] = significance
                    if variant_proportion <= low_level_cutoff:
                        snp['level'] = "low"
                    elif variant_proportion >= high_level_cutoff:
                        snp['level'] = "high"
                if not depth_passed:
                    snp['flag'] = "low coverage"
                snp_list.append(snp)
                #print("Found position of interest %d, reference: %s, distribution:%s" % (position, snp_dict[position][0], base_counter))
            # We've covered it, now remove it from the dict so we can see what we might have missed
            del snp_dict[position]
        elif depth_passed and snp_call and snp_count >= mutdepth and snp_call_proportion >= proportion:
            snp = {'name':'unknown', 'position':str(translated), 'depth':str(column_depth), 'reference':reference_call, 'variant':snp_call, 'basecalls':base_counter}
            if 0 in snp_dict:
                (name, *rest, significance) = snp_dict[0][0]
                snp['name'] = name
                snp['significance'] = significance
            snp_list.append(snp)
            #print("SNP found at position %d: %s->%s" % (position, reference_call, alignment_call))
    
    #Check for any positions_of_interest that weren't covered
    snp_dict.pop(0, None)
    for position in snp_dict.keys():
        for (name, reference, variant, significance) in snp_dict[position]:
            snp = {'name':name, 'position':str(position), 'depth':str(0), 'reference':reference, 'variant':variant}
            snp_list.append(snp)

    if not wholegenome: #If reference is whole genome, none of these are going to make sense, and they will make the output too large
        pileup_dict['consensus_sequence'] = consensus_seq
        pileup_dict['depths'] = ",".join(str(n) for n in depth_array)
        pileup_dict['proportions'] = ",".join(prop_array)
    pileup_dict['breadth'] = str(breadth_positions/amplicon_length * 100)
    pileup_dict['discards'] = ",".join(str(n) for n in discard_array)
    pileup_dict['SNPs'] = snp_list
    pileup_dict['average_depth'] = str(avg_depth_total/avg_depth_positions) if avg_depth_positions else "0"
    return pileup_dict 

def _process_pileup(pileup, amplicon, depth, proportion, mutdepth, offset, wholegenome):
    pileup_dict = {}
    snp_dict = _create_snp_dict(amplicon)
    consensus_seq = ""
    snp_list = []
    breadth_positions = 0
    avg_depth_total = avg_depth_positions = 0
    amplicon_length = len(amplicon.sequence)
    depth_array = ["0"] * amplicon_length
    prop_array = ["0"] * amplicon_length
    for pileupcolumn in pileup:
        depth_array[pileupcolumn.pos] = str(pileupcolumn.n)
        position = pileupcolumn.pos+1
        depth_passed = False
        if pileupcolumn.n > 0: #TODO: This is going to end up being specific to these TB assays, maybe have a clever way to make this line optional
            avg_depth_positions += 1
            avg_depth_total += pileupcolumn.n
        if pileupcolumn.n >= depth:
            breadth_positions += 1
            depth_passed = True
        base_counter = Counter()
        for pileupread in pileupcolumn.pileups:
            if pileupread.is_del:
                base_counter.update("_") # XSLT doesn't like '-' as an attribute name, have to use '_'
            elif pileupread.indel > 0: #This means the next position is an insertion
                #print("Found an insertion at position %d, cigar string: %s" % (position, pileupread.alignment.cigarstring))
                start = pileupread.query_position
                end = pileupread.query_position + pileupread.indel + 1
                base_counter.update({pileupread.alignment.query_sequence[start:end]: 1})
            else:
                base_counter.update(pileupread.alignment.query_sequence[pileupread.query_position])
        #print(base_counter)
        ordered_list = base_counter.most_common()
        alignment_call = ordered_list[0][0]
        alignment_call_proportion = ordered_list[0][1] / pileupcolumn.n
        prop_array[pileupcolumn.pos] = "%.3f" % alignment_call_proportion
        reference_call = amplicon.sequence[pileupcolumn.pos]
        if reference_call == '-':
            reference_call = '_' #Need to use '_' instead of '-' for gaps because of XSLT
        if alignment_call != reference_call:
            snp_call = alignment_call
            snp_count = ordered_list[0][1]
            snp_call_proportion = alignment_call_proportion
        elif len(ordered_list) > 1:
            snp_call = ordered_list[1][0]
            snp_count = ordered_list[1][1]
            snp_call_proportion = ordered_list[1][1] / pileupcolumn.n
        else:
            snp_call = snp_count = snp_call_proportion = None
        consensus_seq += alignment_call if alignment_call_proportion >= proportion else "N"
        translated = offset+position
        if position in snp_dict:
            for (name, reference, variant, significance) in snp_dict[position]:
                snp = {'name':name, 'position':str(translated), 'depth':str(pileupcolumn.n), 'reference':reference, 'variant':variant, 'basecalls':base_counter}
                variant_proportion = base_counter[variant]/pileupcolumn.n
                variant_count = base_counter[variant]
                if variant_proportion >= proportion and variant_count >= mutdepth:
                    snp['significance'] = significance
                    if variant_proportion <= low_level_cutoff:
                        snp['level'] = "low"
                    elif variant_proportion >= high_level_cutoff:
                        snp['level'] = "high"
                if not depth_passed:
                    snp['flag'] = "low coverage"
                snp_list.append(snp)
                #print("Found position of interest %d, reference: %s, distribution:%s" % (position, snp_dict[position][0], base_counter))
            # We've covered it, now remove it from the dict so we can see what we might have missed
            del snp_dict[position]
        elif depth_passed and snp_call and snp_count >= mutdepth and snp_call_proportion >= proportion:
            snp = {'name':'unknown', 'position':str(translated), 'depth':str(pileupcolumn.n), 'reference':reference_call, 'variant':snp_call, 'basecalls':base_counter}
            if 0 in snp_dict:
                (name, *rest, significance) = snp_dict[0][0]
                snp['name'] = name
                snp['significance'] = significance
            snp_list.append(snp)
            #print("SNP found at position %d: %s->%s" % (position, reference_call, alignment_call))

    #Check for any positions_of_interest that weren't covered
    snp_dict.pop(0, None)
    for position in snp_dict.keys():
        for (name, reference, variant, significance) in snp_dict[position]:
            snp = {'name':name, 'position':str(position), 'depth':str(0), 'reference':reference, 'variant':variant}
            snp_list.append(snp)

    if not wholegenome: #If reference is whole genome, none of these are going to make sense, and they will make the output too large
        pileup_dict['consensus_sequence'] = consensus_seq
        pileup_dict['depths'] = ",".join(depth_array)
        pileup_dict['proportions'] = ",".join(prop_array)
    pileup_dict['breadth'] = str(breadth_positions/amplicon_length * 100)
    pileup_dict['SNPs'] = snp_list
    pileup_dict['average_depth'] = str(avg_depth_total/avg_depth_positions) if avg_depth_positions else "0"
    return pileup_dict 

def _write_xml(root, xml_file):
    from xml.dom import minidom
    #logging.debug(ElementTree.dump(root))
    dom = minidom.parseString(ElementTree.tostring(root))
    output = open(xml_file, 'w')
    output.write(dom.toprettyxml(indent="  "))
    output.close()
    return xml_file

def _create_snp_dict(amplicon):
    snp_dict = {}
    for snp in amplicon.SNPs:
        name = snp.name if snp.name else "position of interest"
        if not snp.variant or snp.variant == "any":
            for v in {'A', 'G', 'C', 'T'}:
                if v != snp.reference:
                    if snp.position in snp_dict:
                        snp_dict[snp.position].append((name, snp.reference, v, snp.significance))
                    else:
                        snp_dict[snp.position] = [(name, snp.reference, v, snp.significance)]
        else: 
            if snp.position in snp_dict:
                snp_dict[snp.position].append((name, snp.reference, snp.variant, snp.significance))
            else:
                snp_dict[snp.position] = [(name, snp.reference, snp.variant, snp.significance)]
    return snp_dict

def _add_snp_node(parent, snp):
    snp_attributes = {k:snp[k] for k in ('name', 'position', 'depth', 'reference')}
    snp_node = ElementTree.SubElement(parent, 'snp', snp_attributes)
    base_counter = snp.get('basecalls')
    snpcall = snp['variant']
    depth = int(snp['depth'])
    snpcount = base_counter[snpcall] if base_counter else 0
    percent = snpcount/depth*100 if depth else 0
    snpcall_node = ElementTree.SubElement(snp_node, 'snp_call', {'count':str(snpcount), 'percent':str(percent)})
    snpcall_node.text = snpcall
    if 'significance' in snp or 'flag' in snp:
        significance_node = ElementTree.SubElement(snp_node, 'significance')
        if 'significance' in snp:
            significance_node.text = snp['significance'].message
            if snp['significance'].resistance:
                significance_node.set("resistance", snp['significance'].resistance)
            if 'level' in snp:
                significance_node.set("level", snp['level'])
        if 'flag' in snp:
            significance_node.set('flag', snp['flag'])
    if base_counter:
        ElementTree.SubElement(snp_node, 'base_distribution', {k:str(v) for k,v in base_counter.items()})
    return snp_node

def _process_roi(roi, samdata, amplicon_ref, reverse_comp=False):
    roi_dict = {'region':roi.position_range}
    range_match = re.search('(\d*)-(\d*)', roi.position_range)
    if not range_match:
        return roi_dict
    start = int(range_match.group(1)) - 1
    end = int(range_match.group(2))
    if end < start:
        reverse_comp = True
        start,end = end,start
    expected_length = end - start
    aa_sequence_counter = Counter()
    aa_sequence_counter_temp = Counter()
    nt_sequence_counter = Counter()
    depth = 0
    significant = False
    if not roi.aa_sequence:
        roi.aa_sequence = str(DNA(roi.nt_sequence).translate()).replace('*', 'x')
    for read in samdata.fetch(amplicon_ref, start, end):
        rstart = read.reference_start
        alignment_length = read.get_overlap(start, end)
        #throw out reads that either have gaps in the ROI or don't cover the whole ROI
        if alignment_length != expected_length:
            continue
        if rstart <= start:
            qend = qstart = None
            for (qpos, rpos) in read.get_aligned_pairs():
                if rpos == start:
                    qstart = qpos
                if rpos == end:
                    qend = qpos
            #throw out reads with insertions in the ROI
            if not qend or not qstart or qend-qstart != expected_length:
                continue
            nt_sequence = DNA(read.query_sequence[qstart:qend])
            if reverse_comp:
                nt_sequence = nt_sequence.reverse_complement()
            #scikit-bio doesn't support translating degenerate bases currently, so we will just throw out reads with degenerates for now
            if nt_sequence.has_degenerates(): 
                continue
            aa_sequence = nt_sequence.translate()
            aa_string = str(aa_sequence).replace('*', 'x')
            if aa_string:
                nt_sequence_counter.update([str(nt_sequence)])
                aa_sequence_counter_temp.update([aa_string])
                depth += 1
    if len(aa_sequence_counter_temp) == 0:
        roi_dict['flag'] = "region not found"
        return roi_dict
    else:
        for (aa_string, count) in aa_sequence_counter_temp.most_common():
            num_changes = 0
            for i in range(len(roi.aa_sequence)):
                if len(aa_string) <= i or roi.aa_sequence[i] != aa_string[i]:
                    num_changes += 1 
            aa_sequence_counter[(aa_string, num_changes)] = count

    #This next bit is just being saved for backward compatibility. Should deprecate and remove soon
    (aa_consensus, num_changes) = aa_sequence_counter.most_common(1)[0][0]
    nt_consensus = nt_sequence_counter.most_common(1)[0][0]
    reference = roi.aa_sequence
    consensus = aa_consensus
    if roi.nt_sequence:
        reference = roi.nt_sequence
        consensus = nt_consensus
    roi_dict['most_common_aa_sequence'] = aa_consensus
    roi_dict['most_common_nt_sequence'] = nt_consensus
    roi_dict['reference'] = reference
    roi_dict['changes'] = str(num_changes)
    #End backward compatibility code
    roi_dict['aa_sequence_distribution'] = aa_sequence_counter
    roi_dict['nt_sequence_distribution'] = nt_sequence_counter
    roi_dict['depth'] = str(depth)
    return roi_dict

def _process_roi_SMOR(roi, samdata, amplicon_ref, reverse_comp=False):
    from operator import attrgetter
    roi_dict = {'region':roi.position_range}
    range_match = re.search('(\d*)-(\d*)', roi.position_range)
    if not range_match:
        return roi_dict
    start = int(range_match.group(1)) - 1
    end = int(range_match.group(2))
    if end < start:
        reverse_comp = True
        start,end = end,start
    expected_length = end - start
    aa_sequence_counter = Counter()
    aa_sequence_counter_temp = Counter()
    nt_sequence_counter = Counter()
    depth = 0
    if not roi.aa_sequence:
        roi.aa_sequence = str(DNA(roi.nt_sequence).translate()).replace('*', 'x')
    reads = iter(sorted(samdata.fetch(amplicon_ref, start, end), key=attrgetter('query_name')))
    for read, pair in pairwise(reads):
        if read.query_name != pair.query_name:
            continue
        rstart1 = read.reference_start
        rstart2 = pair.reference_start
        alignment_length1 = read.get_overlap(start, end)
        alignment_length2 = pair.get_overlap(start, end)
        #throw out reads that either have gaps in the ROI or don't cover the whole ROI
        if alignment_length1 != expected_length or alignment_length2 != expected_length:
            continue
        if rstart1 <= start:
            qend = qstart = None
            for (qpos, rpos) in read.get_aligned_pairs():
                if rpos == start:
                    qstart = qpos
                if rpos == end:
                    qend = qpos
            #throw out reads with insertions in the ROI
            if not qend or not qstart or qend-qstart != expected_length:
                continue
            nt_sequence = DNA(read.query_sequence[qstart:qend])
        if rstart2 <= start:
            qend = qstart = None
            for (qpos, rpos) in pair.get_aligned_pairs():
                if rpos == start:
                    qstart = qpos
                if rpos == end:
                    qend = qpos
            #throw out reads with insertions in the ROI
            if not qend or not qstart or qend-qstart != expected_length:
                continue
            nt_sequence2 = DNA(pair.query_sequence[qstart:qend])
        if nt_sequence != nt_sequence2:
            continue
        else:
            if reverse_comp:
                nt_sequence = nt_sequence.reverse_complement()
            #scikit-bio doesn't support translating degenerate bases currently, so we will just throw out reads with degenerates for now
            if nt_sequence.has_degenerates(): 
                continue
            aa_sequence = nt_sequence.translate()
            aa_string = str(aa_sequence).replace('*', 'x')
            if aa_string:
                nt_sequence_counter.update([str(nt_sequence)])
                aa_sequence_counter_temp.update([aa_string])
                depth += 1
    if len(aa_sequence_counter_temp) == 0:
        roi_dict['flag'] = "region not found"
        return roi_dict
    else:
        for (aa_string, count) in aa_sequence_counter_temp.most_common():
            num_changes = 0
            for i in range(len(roi.aa_sequence)):
                if len(aa_string) <= i or roi.aa_sequence[i] != aa_string[i]:
                    num_changes += 1 
            aa_sequence_counter[(aa_string, num_changes)] = count
    
    (aa_consensus, num_changes) = aa_sequence_counter.most_common(1)[0][0]
    nt_consensus = nt_sequence_counter.most_common(1)[0][0]
    reference = roi.aa_sequence
    consensus = aa_consensus
    if roi.nt_sequence:
        reference = roi.nt_sequence
        consensus = nt_consensus
    roi_dict['most_common_aa_sequence'] = aa_consensus
    roi_dict['most_common_nt_sequence'] = nt_consensus
    roi_dict['reference'] = reference
    roi_dict['changes'] = str(num_changes)
    roi_dict['aa_sequence_distribution'] = aa_sequence_counter
    roi_dict['nt_sequence_distribution'] = nt_sequence_counter
    roi_dict['depth'] = str(depth)
    return roi_dict

def _add_roi_node(parent, roi, roi_dict, depth, proportion, mutdepth, smor):
    global low_level_cutoff, high_level_cutoff
    nonsynonymous = False
    if "flag" in roi_dict:
        roi_node = _add_dummy_roi_node(parent, roi)
        significance_node = ElementTree.SubElement(roi_node, "significance", {'flag':roi_dict['flag']})
        if roi.significance.resistance:
            significance_node.set("resistance", roi.significance.resistance)
        return roi_node
    # Smart SMOR -- thresholds and proportion filter are a function of the number of SMOR reads at position
    if smor:
        (proportion, low_level_cutoff, high_level_cutoff) = _compute_thresholds_SMOR(int(roi_dict['depth']))
    roi_attributes = {k:roi_dict[k] for k in ('region', 'reference', 'depth')}
    roi_attributes['name'] = str(roi.name)
    roi_node = ElementTree.SubElement(parent, "region_of_interest", roi_attributes)
    reporting_threshold = max(mutdepth, math.ceil(int(roi_dict['depth']) * proportion))
    dominant_count = 0; #Number of reads containing the most common amino acid sequence
    aa_seq_counter = roi_dict['aa_sequence_distribution']
    for ((seq, aa_changes), count) in aa_seq_counter.most_common():
        if dominant_count == 0:
            dominant_count = count
        if count >= reporting_threshold:
            aa_seq_node = ElementTree.SubElement(roi_node, "amino_acid_sequence", {'count':str(count), 'percent':str(count/int(roi_dict['depth'])*100), 'aa_changes':str(aa_changes)})
            aa_seq_node.text = seq
            if aa_changes > 0:
                nonsynonymous = True
        else:
            break #Since they are returned in order by count, as soon as one is below the threshold the rest will be as well
    nt_seq_counter = roi_dict['nt_sequence_distribution']
    for (seq, count) in nt_seq_counter.most_common():
        if count >= reporting_threshold:
            nt_seq_node = ElementTree.SubElement(roi_node, "nucleotide_sequence", {'count':str(count), 'percent':str(count/int(roi_dict['depth'])*100)})
            nt_seq_node.text = seq
        else:
            break #Since they are returned in order by count, as soon as one is below the threshold the rest will be as well
    low_level = True
    high_level = False
    significant = False
    for mutation in roi.mutations:
        if roi.nt_sequence:
            count = nt_seq_counter[mutation]
            mutant_proportion = count/int(roi_dict['depth'])
        else:
            count = aa_seq_counter[next((k for k in aa_seq_counter.keys() if k[0] == mutation), None)]
            mutant_proportion = count/int(roi_dict['depth'])
        mutation_node = ElementTree.SubElement(roi_node, 'mutation', {'name':str(roi.name)+mutation, 'count':str(count), 'percent':str(mutant_proportion*100)})
        mutation_node.text = mutation
        if mutant_proportion >= proportion and count >= mutdepth:
            significant = True
            if mutant_proportion > low_level_cutoff:
                low_level = False
            if mutant_proportion >= high_level_cutoff:
                high_level = True
    if significant:
        significance_node = ElementTree.SubElement(roi_node, "significance")
        significance_node.text = roi.significance.message
        if roi.significance.resistance:
            significance_node.set("resistance", roi.significance.resistance)
        if int(roi_dict['depth']) < depth:
            significance_node.set("flag", "low coverage")
        if low_level:
            significance_node.set("level", "low")
        elif high_level:
            significance_node.set("level", "high")
    elif len(roi.mutations) == 0 and dominant_count >= mutdepth and (('changes' in roi_dict and int(roi_dict['changes']) > 0) or nonsynonymous):
        significance_node = ElementTree.SubElement(roi_node, "significance", {'changes':roi_dict['changes']})
        significance_node.text = roi.significance.message                       
        if roi.significance.resistance:
            significance_node.set("resistance", roi.significance.resistance)
        if int(roi_dict['depth']) < depth:
            significance_node.set("flag", "low coverage")
    elif int(roi_dict['depth']) < depth: # No significance but still need to flag it for low coverage
        significance_node = ElementTree.SubElement(roi_node, "significance")
        if roi.significance.resistance:
            significance_node.set("resistance", roi.significance.resistance)
        significance_node.set("flag", "low coverage")
    ElementTree.SubElement(roi_node, 'aa_sequence_distribution', {k[0]:str(v) for k,v in aa_seq_counter.items() if v>1}) #key is a tuple of (sequence, changes) and I just want the sequence
    ElementTree.SubElement(roi_node, 'nt_sequence_distribution', {k:str(v) for k,v in nt_seq_counter.items() if v>1})
    return roi_node

def _add_dummy_roi_node(parent, roi):
    reference = roi.aa_sequence
    if roi.nt_sequence:
        reference = roi.nt_sequence
    roi_attributes = {'region':roi.position_range, 'reference':reference, 'depth':"0"}
    roi_node = ElementTree.SubElement(parent, "region_of_interest", roi_attributes)
    for mutation in roi.mutations:
        mutation_node = ElementTree.SubElement(roi_node, 'mutation', {'name':str(roi.name)+mutation, 'count':"0", 'percent':"0"})
        mutation_node.text = mutation
    return roi_node

def _compute_thresholds_SMOR(smor_count):
    proportion = 0
    low_level_cutoff = 0
    high_level_cutoff = 0
    if smor_count >= 5000:
        proportion = 0.001
        low_level_cutoff = 0.01
        high_level_cutoff = 0.5
    elif smor_count >= 500:
        proportion = 0.01
        low_level_cutoff = 0.01
        high_level_cutoff = 0.5
    elif smor_count >= 50:
        proportion = 0.1
        low_level_cutoff = 0.1
        high_level_cutoff = 0.5
    elif smor_count >= 25:
        proportion = 0.2
        low_level_cutoff = 0.2
        high_level_cutoff = 0.5
    return (proportion, low_level_cutoff, high_level_cutoff)

def _merge_reads(read, pair):
    from copy import deepcopy
    rstart = read.query_alignment_start
    rend = read.query_alignment_end
    pstart = pair.query_alignment_start
    pend = pair.query_alignment_end
    merged_read = deepcopy(read)
    if pstart <= rend:  #There is overlap that will need to be processed
        pass
    if pstart > rend+1: #There is a gap that will need to be filled with Ns
        pass
    

def _verify_percent_identity(samdata, ref_name, amplicon, percid, merge):
    temp_file = "%s_%s_temp.bam" % (os.path.splitext(os.path.basename(samdata.filename.decode("utf-8")))[0], ref_name)
    outdata = pysam.AlignmentFile(temp_file, "wb", template=samdata)
    amp_length = len(amplicon.sequence)
    discarded_reads = 0
    seq_counter = Counter()
    aligned_reads = []
    logging.debug("Checking %s for amplicon %s, length %i" % (samdata.filename, ref_name, amp_length))
    if merge:
        reads = iter(sorted(samdata.fetch(ref_name), key=attrgetter('query_name')))
        for read, pair in pairwise(reads):
            if read.alignment.query_name != pair.alignment.query_name:
                continue
            aligned_reads.append(_merge_reads(read, pair))
    else:
        aligned_reads = samdata.fetch(ref_name)
    for read in aligned_reads:
        length = read.infer_query_length(False)
        amp_length = len(amplicon.sequence) #reset amp_length in case we altered it in the last iteration
        logging.debug("Read %s, aligned length %i, total read length %i" % (read.query_name, read.query_alignment_length, length))
        if read.is_unmapped:
            logging.debug("\tRead is unmapped, skipping....");
            continue
        if read.query_alignment_length / length >= percid: #Using length instead of amp_length to compare to query instead of reference
            matches = 0
            for (qpos, rpos, seq) in read.get_aligned_pairs(with_seq=True):
                query = read.query_sequence[qpos] if qpos else "None"
                logging.debug("\tqpos: %i\trpos: %i\tseq: %s\tquery[qpos]: %s" % (qpos or -1, rpos or -1, seq, query)) 
                #if there is a gap in the alignment, extend the length of the query or reference accordingly
                if rpos is None:
                    amp_length += 1
                elif qpos is None:
                    length += 1
                else:
                    if read.query_sequence[qpos].upper() == seq.upper():
                        matches += 1
            if matches / length >= percid: #Using length instead of amp_length to compare to query instead of reference
                logging.debug("\t\tFound %i matches, keeping..." % matches)
                outdata.write(read)
            else:
                logging.debug("\t\tFound %i matches, discarding..." % matches)
                discarded_reads += 1
                seq_counter.update([read.query_sequence])
        else:
            logging.debug("\t\tToo short, discarding...")
            discarded_reads += 1
            seq_counter.update([read.query_sequence])
    outdata.close()
    pysam.index(temp_file)
    return (temp_file, discarded_reads, seq_counter)

class CLIError(Exception):
    '''Generic exception to raise and log different fatal errors.'''
    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg
    def __str__(self):
        return self.msg
    def __unicode__(self):
        return self.msg

def main(argv=None): # IGNORE:C0111
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    if __name__ == '__main__':
        program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    else:
        program_shortdesc = __doc__.split("\n")[1]    
    #program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s

  Created by TGen North on %s.
  Copyright 2015 TGen North. All rights reserved.

  Available for academic and research use only under a license
  from The Translational Genomics Research Institute (TGen)
  that is free for non-commercial use.

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
''' % (program_shortdesc, str(__date__))

    try:
        # Setup argument parser
        parser = argparse.ArgumentParser(description=program_license, formatter_class=argparse.RawDescriptionHelpFormatter)
        required_group = parser.add_argument_group("required arguments")
        required_group.add_argument("-j", "--json", metavar="FILE", required=True, help="JSON file of assay descriptions. [REQUIRED]")
        required_group.add_argument("-b", "--bam", metavar="FILE", required=True, help="BAM file to analyze. [REQUIRED]")
        #required_group.add_argument("-r", "--ref", metavar="FILE", required=True, help="reference fasta file, should already be indexed. [REQUIRED]")
        #parser.add_argument("-o", "--out-dir", dest="odir", metavar="DIR", help="directory to write output files to. [default: `pwd`]")
        required_group.add_argument("-o", "--out", metavar="FILE", required=True, help="XML file to write output to. [REQUIRED]")
        #parser.add_argument("-n", "--name", help="sample name, if not provided it will be derived from BAM file")
        parser.add_argument("-d", "--depth", default=100, type=int, help="minimum read depth required to consider a position covered. [default: 100]")
        parser.add_argument("--breadth", default=0.8, type=float, help="minimum breadth of coverage required to consider an amplicon as present. [default: 0.8]")
        parser.add_argument("-p", "--proportion", default=0.1, type=float, help="minimum proportion required to call a mutation at a given locus. [default: 0.1]")
        parser.add_argument("-m", "--mutation-depth", dest="mutdepth", default=10, type=int, help="minimum number of reads required to call a mutation at a given locus. [default: 10]")
        identity_group = parser.add_argument_group("identity filter options")
        identity_group.add_argument("-i", "--identity", dest="percid", default=0, type=float, help="minimum percent identity required to align a read to a reference amplicon sequence. [default: 0]")
        keep_discarded_group = identity_group.add_mutually_exclusive_group()
        keep_discarded_group.add_argument("-k", "--keep", action="store_true", default=False, help="keep filtered reads. [default: True]")
        keep_discarded_group.add_argument("--no-keep", action="store_false", dest="keep", help="discard filtered reads.")
        merge_group = identity_group.add_mutually_exclusive_group()
        merge_group.add_argument("--merge", action="store_true", default=False, help="merge paired reads. [default: False]")
        merge_group.add_argument("--no-merge", action="store_false", dest="merge", help="do not merge paired reads.")
        parser.add_argument("-s", "--smor", action="store_true", default=False, help="perform SMOR analysis with overlapping reads. [default: False]")
        parser.add_argument("-V", "--version", action="version", version=program_version_message)
        parser.add_argument("-D", "--debug", action="store_true", default=False, help="write <sample_name>.log file with debugging information")
        parser.add_argument("-w", "--whole-genome", action="store_true", dest="wholegenome", default=False, help="JSON file uses a whole genome reference, so don't write out the consensus, depth, and proportion arrays for each sample")
     
        # Process arguments
        args = parser.parse_args()

        json_fp = args.json
        bam_fp = args.bam
        out_fp = args.out
        depth = args.depth
        breadth = args.breadth
        proportion = args.proportion
        mutdepth = args.mutdepth
        percid = args.percid
        keep = args.keep
        merge = args.merge
        smor = args.smor
        debug = args.debug
        wholegenome = args.wholegenome
        #ref_fp = args.ref
        #out_dir = args.odir
        #if not out_dir:
        #    out_dir = os.getcwd()
       
        #out_dir = dispatcher.expandPath(out_dir)
        #if not os.path.exists(out_dir):
        #    os.makedirs(out_dir)

        assay_list = assayInfo.parseJSON(json_fp)
        samdata = pysam.AlignmentFile(bam_fp, "rb")
        #reference = pysam.FastaFile(ref_fp)
        
        sample_dict = {}
        if 'RG' in samdata.header.to_dict() :
            sample_dict['name'] = samdata.header.to_dict()['RG'][0]['ID']
        else:
            sample_dict['name'] = os.path.splitext(os.path.basename(bam_fp))[0]
        sample_dict['mapped_reads'] = str(samdata.mapped)
        sample_dict['unmapped_reads'] = str(samdata.unmapped)
        sample_dict['unassigned_reads'] = str(samdata.nocoordinate)
        sample_dict['depth_filter'] = str(depth)
        #When doing SMOR analysis, proportion filter will be a function of SMOR count at a given postion/ROI,
        # and will ignore passed-in value. Let's specifically set to zero, to make sure all get reported.
        if smor:
            sample_dict['SMOR'] = 'True'
            proportion = 0.001
        sample_dict['proportion_filter'] = str(proportion)
        sample_dict['breadth_filter'] = str(breadth)
        sample_dict['mutation_depth_filter'] = str(mutdepth)
        if percid:
            sample_dict['identity_filter'] = str(percid)
        sample_dict['json_file'] = json_fp
        sample_dict['bam_file'] = bam_fp
        sample_node = ElementTree.Element("sample", sample_dict)

        if debug:
            logfile = "%s.log" % sample_dict['name']

            logging.basicConfig(level=logging.DEBUG,
                                format='%(asctime)s %(levelname)-8s %(message)s',
                                datefmt='%m/%d/%Y %H:%M:%S',
                                filename=logfile,
                                filemode='w')

        #out_fp = os.path.join(out_dir, sample_dict['name']+".xml")
        
        for assay in assay_list:
            assay_dict = {}
            assay_dict['name'] = assay.name
            assay_dict['type'] = assay.assay_type
            assay_dict['function'] = assay.target.function or ""
            assay_dict['gene'] = assay.target.gene_name or ""
            #offset is where the amplicon sits relative to a reference, subtract 1 to make 0-based position
            #report the lesser of start and end in case amplicon is on reverse strand
            try:
                offset = min(int(assay.target.start_position), int(assay.target.end_position))-1
            except:
                offset = 0
            assay_node = ElementTree.SubElement(sample_node, "assay", assay_dict)
            ref_name = assay.name
            reverse_comp = assay.target.reverse_comp
            for amplicon in assay.target.amplicons:
                temp_file = None
                ref_name = assay.name + "_%s" % amplicon.variant_name if amplicon.variant_name else assay.name
                amplicon_dict = {}
                seq_counter = None
                if percid:
                    (temp_file, discarded_reads, seq_counter) = _verify_percent_identity(samdata, ref_name, amplicon, percid, merge)
                    samdata = pysam.AlignmentFile(temp_file, "rb")
                    amplicon_dict['discarded_reads'] = str(discarded_reads)
                elif samdata.closed:
                    samdata = pysam.AlignmentFile(bam_fp, "rb")
                amplicon_dict['reads'] = str(samdata.count(ref_name))
                if amplicon.variant_name:
                    amplicon_dict['variant'] = amplicon.variant_name
                amplicon_node = ElementTree.SubElement(assay_node, "amplicon", amplicon_dict)
                if keep and seq_counter:
                    ElementTree.SubElement(amplicon_node, "sequence_distribution", {k:str(v) for k,v in seq_counter.items()})
                if samdata.count(ref_name) == 0:
                    significance_node = ElementTree.SubElement(amplicon_node, "significance", {"flag":"no coverage"})
                    #Check for indeterminate resistances
                    resistances = set()
                    if amplicon.significance and amplicon.significance.resistance:
                        resistances.add(amplicon.significance.resistance)
                    for snp in amplicon.SNPs:
                        name = snp.name if snp.name else "position of interest"
                        dummy_snp = {'name':name, 'position':str(snp.position), 'depth':"0", 'reference':snp.reference, 'variant':snp.variant, 'basecalls':None}
                        _add_snp_node(amplicon_node, dummy_snp)
                        if snp.significance.resistance:
                            resistances.add(snp.significance.resistance)
                    for roi in amplicon.ROIs:
                        _add_dummy_roi_node(amplicon_node, roi)
                        if roi.significance.resistance:
                            resistances.add(roi.significance.resistance)
                    if resistances:        
                        significance_node.set("resistance", ",".join(resistances))
                else:
                    if amplicon.significance or samdata.count(ref_name) < depth:
                        significance_node = ElementTree.SubElement(amplicon_node, "significance")
                        if amplicon.significance:
                            significance_node.text = amplicon.significance.message
                            if amplicon.significance.resistance:
                                significance_node.set("resistance", amplicon.significance.resistance)
                        if samdata.count(ref_name) < depth:
                            significance_node.set("flag", "low coverage")
                            #Check for indeterminate resistances
                            resistances = set()
                            if amplicon.significance and amplicon.significance.resistance:
                                resistances.add(amplicon.significance.resistance)
                            for snp in amplicon.SNPs:
                                if snp.significance.resistance:
                                    resistances.add(snp.significance.resistance)
                            for roi in amplicon.ROIs:
                                if roi.significance.resistance:
                                    resistances.add(roi.significance.resistance)
                            if resistances:        
                                significance_node.set("resistance", ",".join(resistances))

                    pileup = samdata.pileup(ref_name, max_depth=1000000)
                    if smor:
                        amplicon_data = _process_pileup_SMOR(pileup, amplicon, depth, proportion, mutdepth, offset, wholegenome)
                    else:
                        amplicon_data = _process_pileup(pileup, amplicon, depth, proportion, mutdepth, offset, wholegenome)
                    if float(amplicon_data['breadth']) < breadth*100:
                        significance_node = amplicon_node.find("significance")
                        if significance_node is None:
                            significance_node = ElementTree.SubElement(amplicon_node, "significance")
                        if not significance_node.get("flag"):
                            significance_node.set("flag", "insufficient breadth of coverage")
                    for snp in amplicon_data['SNPs']:
                        _add_snp_node(amplicon_node, snp)
                        # This would be helpful, but count_coverage is broken in python3
                        #print(samdata.count_coverage(ref_name, snp.position-1, snp.position))
                    del amplicon_data['SNPs']
                    _write_parameters(amplicon_node, amplicon_data)

                    for roi in amplicon.ROIs:
                        if smor:
                            roi_dict = _process_roi_SMOR(roi, samdata, ref_name, reverse_comp)
                        else:
                            roi_dict = _process_roi(roi, samdata, ref_name, reverse_comp)
                        _add_roi_node(amplicon_node, roi, roi_dict, depth, proportion, mutdepth, smor)

                if temp_file:
                    samdata.close()
                    os.remove(temp_file)
                    os.remove(temp_file+".bai")
                    samdata = pysam.AlignmentFile(bam_fp, "rb")

        if samdata.is_open():
            samdata.close()
        _write_xml(sample_node, out_fp)

        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception as e:
        if DEBUG or TESTRUN:
            raise(e)
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2

if __name__ == "__main__":
    if DEBUG:
        pass
    if TESTRUN:
        import doctest
        doctest.testmod()
    if PROFILE:
        import cProfile
        import pstats
        profile_filename = 'asap.analyzeAmplicons_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())

#!/usr/bin/env python3
# encoding: utf-8
'''
asap.bamProcessor -- Process BAM alignment files with an AssayInfo JSON file and generate XML for the results

asap.bamProcessor

@author:     Darrin Lemmer

@copyright:  2015,2025 TGen North. All rights reserved.

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
from skbio import TabularMSA, DNA
from skbio.alignment import local_pairwise_align_ssw

from asap import assayInfo
from asap import __version__
from asap import cmdParser
# https://github.com/martinblech/xmltodict
import json
import xmltodict
import numpy as np
import array as arr


__all__ = []
__updated__ = '2025-05-16'
__date__ = '2015-07-16'

DEBUG = 1
INFO = 0
TESTRUN = 0
PROFILE = 0
REMOVE_TEMP = True

low_level_cutoff = 0.01
high_level_cutoff = 0.50
proportion = 0.1

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

def _process_pileup(pileup, amplicon, depth, proportion, mutdepth, offset, wholegenome, base_qual, con_prop, fill_gap_char, fill_del_char):
    global low_level_cutoff, high_level_cutoff
    pileup_dict = {}
    snp_dict = _create_snp_dict(amplicon)
    deletion_counter = Counter() #keep track of deletions by read name
    consensus_seq = ""
    if fill_gap_char:
        gapfilled_consensus_seq = ""
    snp_list = []
    breadth_positions = 0
    avg_depth_total = avg_depth_positions = 0
    amplicon_length = len(amplicon.sequence)
    depth_array = [0] * amplicon_length
    quality_discard_array = [0] * amplicon_length
    prop_array = ["0"] * amplicon_length
    n_read_array = [0] * amplicon_length # New array to count 'N' reads
    previous_position = 0
    # for each position in alignment/pileup
    for pileupcolumn in pileup:
        base_counter = Counter()
        position = pileupcolumn.pos+1
        # This fills gaps in the alignment with n's or user defined char
        if fill_gap_char:
            if previous_position+1 < position: #We've skipped some positions in the alignment
                #print("%i, %i" % (previous_position, position))
                for i in range(previous_position+1, position):
                    gapfilled_consensus_seq += fill_gap_char #Fill in the gap
        previous_position = position
        depth_array[pileupcolumn.pos] = pileupcolumn.n
        depth_passed = False
        passed_Qual_filter = 0
        for pileupread in pileupcolumn.pileups:
            #print("processing read, qual=%i" % pileupread.alignment.query_qualities[pileupread.query_position])
            try:
                # Tanner: Check for 'N' bases first
                if pileupread.query_position is not None and pileupread.alignment.query_sequence[pileupread.query_position].upper() == 'N':
                    n_read_array[pileupcolumn.pos] += 1
                    continue
                if pileupread.is_del:
                    #This position in the alignment is a deletion in the query sequence, therefore it has no quality score
                    # Let's use the average of the quality scores of the two aligned bases flanking the deletion
                    qscore = (pileupread.alignment.query_qualities[pileupread.query_position_or_next] +
                              pileupread.alignment.query_qualities[pileupread.query_position_or_next - 1]) / 2
                    if qscore >= base_qual:
                        passed_Qual_filter += 1
                        base_counter.update({"_" : 1})
                    else:
                        quality_discard_array[pileupcolumn.pos] += 1 # TP Updated: previous quality_discard_array[pileupcolumn.pos] 
                elif pileupread.alignment.query_qualities[pileupread.query_position] >= base_qual: # check here
                    passed_Qual_filter += 1
                    if pileupread.indel < 0: #This means the next position is a deletion, we'll process later
                        for d in range(1, abs(pileupread.indel)+1):
                            deletion_counter.update({str(position + d)})
                    if pileupread.indel > 0: #This means the next position is an insertion, unlike with deletions, this we can process now
                        start = pileupread.query_position
                        end = pileupread.query_position + pileupread.indel + 1
                        base_counter.update({pileupread.alignment.query_sequence[start:end]: 1})
                    else:
                        base_counter.update(pileupread.alignment.query_sequence[pileupread.query_position])
                else:
                    quality_discard_array[pileupcolumn.pos] += 1
            except Exception as e:
                if str(e.__class__.__name__) != "TypeError":
                    print("Unexpected error:", sys.exc_info()[0])
                    pass
                quality_discard_array[pileupcolumn.pos] += 1 #check here
                pass

        column_depth = passed_Qual_filter #check this, this will count bases that have been filtered out by quality?
        depth_array[pileupcolumn.pos] = passed_Qual_filter #reset to depth that passed qual filter

        if column_depth > 0: #TODO: This is going to end up being specific to these TB assays (with flanking sequence), maybe have a clever way to make this line optional
            avg_depth_positions += 1
            avg_depth_total += column_depth
        if column_depth >= depth:
            breadth_positions += 1
            depth_passed = True
        ordered_list = base_counter.most_common()
        if not ordered_list: #No coverage, should only happen here if all reads were thrown out because of quality
            consensus_seq += "N"
            if fill_gap_char:
                gapfilled_consensus_seq += "N"
            continue
        alignment_call = ordered_list[0][0]
        alignment_call_proportion = ordered_list[0][1] / column_depth
        prop_array[pileupcolumn.pos] = "%.3f" % alignment_call_proportion
        reference_call = amplicon.sequence[pileupcolumn.pos]
        if reference_call == '-':
            reference_call = '_' #Need to use '_' instead of '-' for gaps because of XSLT
        #if alignment_call != reference_call:
        #    snp_call = alignment_call
        #    snp_count = ordered_list[0][1]
        #    snp_call_proportion = alignment_call_proportion
        #elif len(ordered_list) > 1:
        #    snp_call = ordered_list[1][0]
        #    snp_count = ordered_list[1][1]
        #    snp_call_proportion = ordered_list[1][1] / column_depth
        # Initialize SNP variables TP added 2025 ########
        snp_call = None
        snp_count = None
        snp_call_proportion = None
        # Find the first valid SNP candidate (a non-reference, non-ambiguous base)
        for base, count in ordered_list:
            if base != reference_call and base.upper() in {'A', 'C', 'G', 'T'}:
                snp_call = base
                snp_count = count
                snp_call_proportion = count / column_depth
                break # Exit the loop once a valid SNP is found
        else:
        #    snp_call = snp_count = snp_call_proportion = None
        # This 'else' block executes if the loop completes without finding a valid SNP
            snp_call = snp_count = snp_call_proportion = None
        #Generate consensus call at this pos
        #consensus_seq += alignment_call if alignment_call_proportion >= consensus_proportion else "N"
        # unless the alignment_call is a deletion, and >50% -- don't ever replace deletions with Ns
        # or if coverage is less than the depth threshold, then always call N
        if not depth_passed: # N's if we don't have enough coverage
            consensus_seq += "N"
            if fill_gap_char:
                gapfilled_consensus_seq += "N"
        elif alignment_call != "_":
            if alignment_call_proportion >= con_prop:
                consensus_seq += alignment_call
                if fill_gap_char:
                    gapfilled_consensus_seq += alignment_call
            else: #Consensus proportion not high enough
                consensus_seq += "N"
                if fill_gap_char:
                    gapfilled_consensus_seq += "N"
        else:
            if alignment_call_proportion <= 0.5: #Verify that the gap call is truly greater than 50%
                consensus_seq += "N"
                if fill_gap_char:
                    gapfilled_consensus_seq += "N"
            else:
                if fill_del_char: #Put in gaps if user requested them
                    consensus_seq += fill_del_char
                    if fill_gap_char:
                        gapfilled_consensus_seq += fill_del_char

        if position >= abs(offset) and offset < 0: #if the offset is negative, ie. amplicon starts before beginning of the gene, then when converting to gene-based coordinates need to make offset 1 unit more positive to account for there being no 0-base in gene-coordinates
            translated = position + (offset + 1)
        else:
            translated = position + offset #normal case where gene encompasses the amplicon
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
            # We've covered it, now remove it from the dict so we can see what we might have missed
            del snp_dict[position]
        elif depth_passed and snp_call and snp_count >= mutdepth and snp_call_proportion >= proportion:
            snp = {'name':'unknown', 'position':str(translated), 'depth':str(column_depth), 'reference':reference_call, 'variant':snp_call, 'basecalls':base_counter}
            if 0 in snp_dict:
                (name, *rest, significance) = snp_dict[0][0]
                snp['name'] = name
                snp['significance'] = significance
            snp_list.append(snp)
    #Check for any positions_of_interest that weren't covered
    snp_dict.pop(0, None)
    for position in snp_dict.keys():
        for (name, reference, variant, significance) in snp_dict[position]:
            snp = {'name':name, 'position':str(position), 'depth':str(0), 'reference':reference, 'variant':variant}
            snp_list.append(snp)
    if not wholegenome: #If reference is whole genome, none of these are going to make sense, and they will make the output too large
        pileup_dict['consensus_sequence'] = consensus_seq
        if fill_gap_char:
            pileup_dict['gapfilled_consensus_sequence'] = gapfilled_consensus_seq
        pileup_dict['depths'] = ",".join(str(n) for n in depth_array)
        pileup_dict['proportions'] = ",".join(prop_array)
        pileup_dict['n_reads'] = ",".join(str(n) for n in n_read_array) # New: Add the 'N' read counts
    pileup_dict['breadth'] = str(breadth_positions/amplicon_length * 100)
    pileup_dict['quality_discards'] = ",".join(str(n) for n in quality_discard_array)
    pileup_dict['SNPs'] = snp_list
    pileup_dict['average_depth'] = str(avg_depth_total/avg_depth_positions) if avg_depth_positions else "0"
    return pileup_dict

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

def _process_roi(roi, samdata, amplicon_ref, amplicon_ref_len, reverse_comp=False):
    from operator import attrgetter
    roi_dict = {'region':roi.position_range}
    range_list = roi.position_range.split(",")
    roi_dict['errors'] = {}
    aa_sequence_counter = Counter()
    aa_sequence_counter_temp = Counter()
    nt_sequence_counter = Counter()
    depth = 0
    for pos_range in range_list:
        print("Handling position: ", pos_range)
        range_match = re.search('(\d*)-(\d*)', pos_range)
        if not range_match:
            continue
        start = int(range_match.group(1)) - 1
        end = int(range_match.group(2))
        if end < start:
            reverse_comp = True
            start,end = end,start
        expected_length = end - start
        #check if the roi spans the whole reference, if so can use .query_alignment_sequence to get the whole sequence without running into the problem of the loop not getting to last base
        #still will have a problem when start != 0 but end == amplicon_ref_len
        use_query_alignment_seq = False
        if end == amplicon_ref_len and start == 0:
            use_query_alignment_seq = True
        aligned_reads = sorted(samdata.fetch(amplicon_ref, start, end), key=attrgetter('query_name'))
        big_reads = []
        #check if reads are long enough, if not then merge
        n = failed = 0
        for read, pair in pairwise(aligned_reads):
            if read.reference_end == None or read.reference_start == None:
                continue
            n += 1
            if read.reference_end - read.reference_start < (expected_length * 0.9):
                failed +=1
        if n == 0:
            # roi_dict['flag'] = "region not found"
            print("WARNING::region not found")
            roi_dict['errors'][pos_range] = "region not found"
            continue
            # return roi_dict
        else: # let's merge
            proportion_failed = failed/n
            if proportion_failed >= .95:
                logging.debug("reads are not as big as roi, merging...")
                reads = sorted(samdata.fetch(amplicon_ref, start, end), key=attrgetter('query_name'))
                big_reads = _process_merge(reads, start ,end)

        if not roi.aa_sequence:
            roi.aa_sequence = str(DNA(roi.nt_sequence).translate()).replace('*', 'x')
        if big_reads != []: #merged
            for read in big_reads:
                nt_sequence = DNA(read)
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
        else:
            aligned_reads = samdata.fetch(amplicon_ref, start, end)
            for read in aligned_reads:
                rstart = read.reference_start
                rend = read.reference_end
                #alignment_length = read.get_overlap(start, end)
                #throw out reads that either have gaps in the ROI or don't cover the whole ROI
                #if alignment_length != expected_length:
                #    continue
                #Keep reads with indels, but do throw out reads that don't cover whole ROI
                if not rend or rstart > start or rend < end:
                    continue
                if rstart <= start:
                    if not use_query_alignment_seq:
                        qend = qstart = None
                        for (qpos, rpos) in read.get_aligned_pairs():
                            if rpos == start:
                                qstart = qpos
                            if rpos == end:
                                qend = qpos
                        nt_sequence = DNA(read.query_sequence[qstart:qend])
                    else:
                        #the ROI is the whole ref so can use .query_alignment_sequence
                        nt_sequence = DNA(read.query_alignment_sequence)
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
        pass #End of loop over ranges
    #Fill out dictionary
    if len(aa_sequence_counter_temp) == 0:
        roi_dict['flag'] = "no regions found"
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

def _add_roi_node(parent, roi, roi_dict, depth, proportion, mutdepth, offset, allele_min_reads):
    global low_level_cutoff, high_level_cutoff
    nonsynonymous = False
    if "flag" in roi_dict:
        roi_node = _add_dummy_roi_node(parent, roi)
        significance_node = ElementTree.SubElement(roi_node, "significance", {'flag':roi_dict['flag']})
        if roi.significance.resistance:
            significance_node.set("resistance", roi.significance.resistance)
        return roi_node
    roi_attributes = {k:roi_dict[k] for k in ('region', 'reference', 'depth')}
    roi_attributes['name'] = str(roi.name)
    roi_node = ElementTree.SubElement(parent, "region_of_interest", roi_attributes)
    if roi_dict["errors"] != {}:
        for postion in roi_dict["errors"].keys():
            err_node = ElementTree.SubElement(roi_node, "error")
            err_node.set("position", position)
            err_node.set("message", roi_dict["errors"][postion])
            pass
        pass
    if not roi.aa_sequence:
        roi.aa_sequence = str(DNA(roi.nt_sequence).translate()).replace('*', 'x')
    roi_node.set('aa_reference', roi.aa_sequence)
    reporting_threshold = max(mutdepth, math.ceil(int(roi_dict['depth']) * proportion))
    #print(proportion, low_level_cutoff, high_level_cutoff, int(roi_dict['depth']), reporting_threshold)
    cutOff = int(roi_dict['depth']) * .02
    range_match = re.search('(\d*)-(\d*)', roi.position_range)
    start = int(range_match.group(1)) - 1
    end = int(range_match.group(2))
    if end < start:
        reverse_comp = True
        start,end = end,start
    dominant_count = 0; #Number of reads containing the most common amino acid sequence
    aa_seq_counter = roi_dict['aa_sequence_distribution']
    aa_allele_count = 0
    #calculate offsets depending on if in positive region of gene or negative
    #adding one when in negative to keep consistent with _process_pileup
    aa_offset_pos = math.floor(offset/3)
    aa_offset_neg = math.floor((offset+1)/3)
    for ((seq, aa_changes), count) in aa_seq_counter.most_common():
        if dominant_count == 0:
            dominant_count = count
        if count >= reporting_threshold:
            aa_seq_node = ElementTree.SubElement(roi_node, "amino_acid_sequence", {'count':str(count), 'percent':str(count/int(roi_dict['depth'])*100), 'aa_changes':str(aa_changes)})
            aa_seq_node.text = seq
            if aa_changes > 0:
                nonsynonymous = True
            #get string of the aa changes
            changes = []
            all_changes = []
            for i in range(len(roi.aa_sequence)):
                if i > len(seq) - 1:
                    change = [i, roi.aa_sequence[i], "_"]
                    all_changes.append(change)
                    changes.append(change)
                elif roi.aa_sequence[i] != seq[i]:
                    change = [i, roi.aa_sequence[i], seq[i]]
                    all_changes.append(change)
                    changes.append(change)
            if len(seq) > len(roi.aa_sequence):
                for i in range(len(roi.aa_sequence), len(seq)):
                    change = [i, "_", seq[i]]
                    all_changes.append(change)
                    changes.append(change)
            #check to see if aa changes are a result of an indel, and if so remove them
            start_of_run = _sequential(changes, 0)
            changes = changes[0:start_of_run]
            shift = 0
            #create change string with '1' and '2' that will be replaced by <b><u> and </u></b> in post-processing
            for change in all_changes:
                loc = change[0] + shift
                #check if change is past last base in seq => an indel at the end of seq
                if loc >= len(seq):
                    temp = seq + '1' + '_' + '2'
                else:
                    temp = seq[0:loc] + '1' + seq[loc] + '2' + seq[loc + 1:]
                shift += 2
                seq = temp
            aa_seq_node.set('underline_seq', seq)
            #create changes strings, adjusting the aa coordinates to be gene-relative
            if all_changes != []:
                change_string = ""
                for all_change in all_changes:
                    if all_change[0] >= abs(aa_offset_pos) and aa_offset_pos < 0: #if the offset is negative, ie. amplicon starts before beginning of the gene, then when converting to gene-based coordinates need to make offset 1 unit more positive to account for there being no 0-base in gene-coordinates
                        all_change[0] = all_change[0] + aa_offset_neg
                    else:
                        all_change[0] = all_change[0] + aa_offset_pos #normal case where gene encompasses the amplicon
                    change_string += ', ' + all_change[1] + str(all_change[0]) + all_change[2]
                aa_seq_node.set('aa_changes_specific_all', change_string)
                change_string = ""
                for change in changes:
                    #don't need to shift the aa coordinates here again because all_changes and changes are filled with same, shallow copied, lists
                    if change[0] < 0:
                        continue
                    else:
                        change_string += change[1] + str(change[0]) + change[2] + ', '
                aa_seq_node.set('aa_changes_specific', change_string)
        else:
            break #Since they are returned in order by count, as soon as one is below the threshold the rest will be as well

    nt_seq_counter = roi_dict['nt_sequence_distribution']
    skbio_reference = None
    if roi.nt_sequence:
        skbio_reference = DNA(roi.nt_sequence)
    for (seq, count) in nt_seq_counter.most_common():
        if count >= reporting_threshold:
            nt_seq_node = ElementTree.SubElement(roi_node, "nucleotide_sequence", {'count':str(count), 'percent':str(count/int(roi_dict['depth'])*100)})
            nt_seq_node.text = seq
            if skbio_reference:
                #align to reference
                alignment, score, start_end_positions = local_pairwise_align_ssw(skbio_reference,DNA(seq))
                #get string of the nt changes
                changes = []
                all_changes = []
                i=1
                for aligned_seq_pos in alignment.iter_positions():
                    if str(aligned_seq_pos[1]) == '-' :
                        change = [i+start, str(aligned_seq_pos[0]), "_"]
                        all_changes.append(change)
                        changes.append(change)
                        i = i+1
                    elif str(aligned_seq_pos[0]) == '-':
                        change = [i+start, "_", str(aligned_seq_pos[1])]
                        all_changes.append(change)
                        changes.append(change)
                        #Don't advance the counter for gaps in the reference
                    elif aligned_seq_pos[0] != aligned_seq_pos[1]:
                        change = [i+start, str(aligned_seq_pos[0]), str(aligned_seq_pos[1])]
                        all_changes.append(change)
                        changes.append(change)
                        i = i+1
                    else:
                        i = i+1
                shift = 0
                #create change string with '1' and '2' that will be replaced by <b><u> and </u></b> in post-processing
                for change in all_changes:
                    loc = change[0] - start - 1 + shift
                    temp = seq
                    if change[2] == '_':
                        temp = seq[0:loc] + '1' + '_' + '2' + seq[loc:]
                    elif change[1] == '_':
                        #What do we do with insertions
                        continue
                    else:
                        temp = seq[0:loc] + '1' + seq[loc] + '2' + seq[loc + 1:]
                    shift += 2
                    seq = temp
                nt_seq_node.set('underline_seq', seq)
                #create changes strings, adjusting the nt coordinates to be gene-relative
                if all_changes != []:
                    change_string = ""
                    for all_change in all_changes:
                        if all_change[0] >= abs(offset) and offset < 0: #if the offset is negative, ie. amplicon starts before beginning of the gene, then when converting to gene-based coordinates need to make offset 1 unit more positive to account for there being no 0-base in gene-coordinates
                            all_change[0] = all_change[0] + (offset + 1)
                        else:
                            all_change[0] = all_change[0] + offset #normal case where gene encompasses the amplicon
                        change_string += ', ' + all_change[1] + str(all_change[0]) + all_change[2]
                    nt_seq_node.set('nt_changes_specific_all', change_string)
                    change_string = ""
                    for change in changes:
                        #don't need to shift the nt coordinates here again because all_changes and changes are filled with same, shallow copied, lists
                        if change[0] < 0:
                            continue
                        else:
                            change_string += change[1] + str(change[0]) + change[2] + ', '
                    nt_seq_node.set('nt_changes_specific', change_string)
            else:
                continue #If we don't have a reference NT sequence, then we can't display the changes
        else:
            break #Since they are returned in order by count, as soon as one is below the threshold the rest will be as well

    allele_count = 0
    for (seq, count) in nt_seq_counter.most_common():
        #get most frequent alleles that have a freq of > 2% (this is an arbitrary cut-off)
        if count >= cutOff:
            allele_node = ElementTree.SubElement(roi_node, "allele_sequence", {'count':str(count), 'percent':str(count/int(roi_dict['depth'])*100),'hash':str(hash(seq))})
            allele_node.text = seq
            allele_count += 1
        else:
            if allele_count < 2:
                allele_count += 1
                allele_node = ElementTree.SubElement(roi_node, "allele_sequence", {'count':str(count), 'percent':str(count/int(roi_dict['depth'])*100),'hash':str(hash(seq))})
                allele_node.text = seq
            else:
                break
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
        if int(roi_dict['depth']) < depth:# that is so incredibly strange that i am typing that habitualy
            significance_node.set("flag", "low coverage")
    elif int(roi_dict['depth']) < depth: # No significance but still need to flag it for low coverage
        significance_node = ElementTree.SubElement(roi_node, "significance")
        if roi.significance.resistance:
            significance_node.set("resistance", roi.significance.resistance)
        significance_node.set("flag", "low coverage")
    #keep all alleles until this point so proportional calculations are correct
    #do not output alleles that have less than allele_min_reads # of reads
    ElementTree.SubElement(roi_node, 'aa_sequence_distribution', {k[0]:str(v) for k,v in aa_seq_counter.items() if v >= allele_min_reads}) #key is a tuple of (sequence, changes) and I just want the sequence
    ElementTree.SubElement(roi_node, 'nt_sequence_distribution', {k:str(v) for k,v in nt_seq_counter.items() if v >= allele_min_reads})
    return roi_node

#returns the index of the start of a run of sequential numbers from some point in array to the end of array
#if such a run does not exist it returns the index of last element + 1
#allows for one element gaps, ie. [1,3,4,5,6] would return 0 but [1,4,5,6] would return 1
def _sequential(arr, n):
    if n == len(arr) - 1:
        return n + 1
    for i in range(n, len(arr) - 1):
        if arr[i][0] + 1 != arr[i+1][0] and arr[i][0] + 2 != arr[i+1][0]:
            return _sequential(arr, i + 1)
    return n

def _add_dummy_roi_node(parent, roi):
    reference = roi.aa_sequence
    if roi.nt_sequence:
        reference = roi.nt_sequence
    roi_attributes = {'region':roi.position_range, 'name':str(roi.name), 'reference':reference, 'depth':"0"}
    roi_node = ElementTree.SubElement(parent, "region_of_interest", roi_attributes)
    for mutation in roi.mutations:
        mutation_node = ElementTree.SubElement(roi_node, 'mutation', {'name':str(roi.name)+mutation, 'count':"0", 'percent':"0"})
        mutation_node.text = mutation
    return roi_node

def _process_merge(reads, start, end):
    big_aligned_reads = []
    for read, pair in pairwise(reads):
        if read.query_name != pair.query_name:
            continue
        refstart1 = read.reference_start
        refend1 = read.reference_end
        refstart2 = pair.reference_start
        refend2 = pair.reference_end
        if refstart1 == None or refend1 == None or refstart2 == None or refend2 == None:
            continue
        #get the farthest left and right positions that either read align to reference
        refstart = refstart1 if refstart1 < refstart2 else refstart2
        refend = refend2 if refend2 > refend1 else refend1
        combined_read = ""
        for ref_pos in range(refstart, refend):
            read_base = pair_base = None
            for (read_qpos, read_rpos) in read.get_aligned_pairs():
                if read_rpos == ref_pos and read_qpos != None:
                    read_base = str(DNA(read.query_sequence[int(read_qpos)]))
                    read_qual = read.query_qualities[int(read_qpos)]
                    break
            for (pair_qpos, pair_rpos) in pair.get_aligned_pairs():
                if pair_rpos == ref_pos and pair_qpos != None:
                    pair_base = str(DNA(pair.query_sequence[int(pair_qpos)]))
                    pair_qual = pair.query_qualities[int(pair_qpos)]
                    break
            if read_base == None and pair_base == None:
                combined_read += 'N'
            elif read_base != None and pair_base == None:
                combined_read += read_base
            elif read_base == None and pair_base != None:
                combined_read += pair_base
            else:
                if read_qual >= pair_qual:
                    combined_read += read_base
                else:
                    combined_read += pair_base
        if len(combined_read) == end-start:
            big_aligned_reads.append(combined_read)
    return big_aligned_reads

class CLIError(Exception):
    '''Generic exception to raise and log different fatal errors.'''
    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg
    def __str__(self):
        return self.msg
    def __unicode__(self):
        return self.msg


def evaluateItem(item, sample_node):
    foundNode = findItem(item, sample_node)
    if not foundNode:
        print("Could not find node for the item in sample node, valued as False")
        return False
        pass
    print("Found Node: ", foundNode)

    if item.evaluation == "depth_greater_than":
        if "depth" in foundNode.attrib:
            return float(foundNode.attrib["depth"]) > float(item.value)
            pass
        else:
            return False
        pass
    if item.evaluation == "depth_less_than":
        if "depth" in foundNode.attrib:
            return float(foundNode.attrib["depth"]) < float(item.value)
            pass
        else:
            return False
        pass

    if item.item_type == "snp":
        if "percentage" in item.evaluation or "count" in item.evaluation:
            for call in foundNode:
                if call.tag == "snp_call":
                    foundNode = call
                    print("Looking at snp_call: ", foundNode)
                    break # This should ensure the topmost snp_call
                    pass

    if item.evaluation == "percentage_greater_than":
        if "percentage" in foundNode.attrib:
            return float(foundNode.attrib["percentage"]) > float(item.value)
            pass
        else:
            return False
        pass
    if item.evaluation == "percentage_less_than":
        if "percentage" in foundNode.attrib:
            return float(foundNode.attrib["percentage"]) < float(item.value)
            pass
        else:
            return False
        pass

    if item.evaluation == "count_greater_than":
        if "count" in foundNode.attrib:
            return float(foundNode.attrib["count"]) > float(item.value)
            pass
        else:
            return False
        pass
    if item.evaluation == "count_less_than":
        if "count" in foundNode.attrib:
            return float(foundNode.attrib["count"]) < float(item.value)
            pass
        else:
            return False
        pass
    # Evaluation empty is simply true
    if not item.evaluation or item.evaluation == "" or str(item.evaluation).lower() == "none":
        return True
        pass
    print("Item contains invalid evaluation: ", item.evaluation)
    return False

def findItem(item, sample_node):
    foundNode = None
    for childEle in sample_node:
        if childEle.tag == "assay":
            if item.item_type.lower() == "assay":
                if item.identity_key not in childEle.attrib:
                    continue
                    pass
                if str(childEle.attrib[item.identity_key]).lower() == item.identity_value.lower():
                    return childEle
                    pass
                pass
        for possibleAmpEle in childEle:
            if possibleAmpEle.tag == "amplicon":
                if item.item_type.lower() == "amplicon":
                    if item.identity_key not in possibleAmpEle.attrib:
                        continue
                        pass
                    if str(possibleAmpEle.attrib[item.identity_key]).lower() == item.identity_value.lower():
                        return possibleAmpEle
                        pass
                    pass
                for possibleSubEle in possibleAmpEle:
                    if possibleSubEle.tag == "roi":
                        if item.item_type.lower() == "roi":
                            if item.identity_key not in possibleSubEle.attrib:
                                continue
                                pass
                            if str(possibleSubEle.attrib[item.identity_key]).lower() == item.identity_value.lower():
                                foundNode = possibleSubEle
                                pass
                            pass
                        if item.item_type.lower() == "mutation":
                            for possibleMutEle in possibleSubEle:
                                if possibleMutEle.tag == "mutation":
                                    if item.identity_key == "base" or item.identity_key == "bases" or item.identity_key == "text" or item.identity_key == "mutation":
                                        if str(possibleMutEle.text).lower() == item.identity_value.lower():
                                            foundNode = possibleMutEle
                                            pass
                                        pass
                                    if item.identity_key not in possibleMutEle.attrib:
                                        continue
                                        pass
                                    if str(possibleMutEle.attrib[item.identity_key]).lower() == item.identity_value.lower():
                                        foundNode = possibleMutEle
                                        pass
                                    pass
                                pass
                            pass
                        pass
                    if possibleSubEle.tag == "snp":
                        if item.item_type.lower() == "snp":
                            if item.identity_key not in possibleSubEle.attrib:
                                continue
                                pass
                            if str(possibleSubEle.attrib[item.identity_key]).lower() == item.identity_value.lower():
                                foundNode = possibleSubEle
                                pass
                            pass
                        pass
                    pass
                pass
            pass
        pass
    #End of findItem func
    return foundNode

def evaluateOperation(operation, sample_node):
    print("Evaluating Operation: ", operation)
    print("On sample: ", sample_node)
    boolList = []
    for child in operation.children:
        if isinstance(child, assayInfo.ITEM):
            print("Evaluating Truthness of Item: ", child)
            boolList.append(evaluateItem(child, sample_node))
            pass
        if isinstance(child, assayInfo.Operation):
            print("Nested Operation: ", child)
            boolList.append(evaluateOperation(child, sample_node))
            pass
        pass
    print("Evaled Bools: ", boolList)
    if operation.operation_type == "AND":
        endBool = True
        for evaldBool in boolList:
            endBool = evaldBool and endBool
            pass
        return endBool
        pass
    if operation.operation_type == "OR":
        endBool = False
        for evaldBool in boolList:
            endBool = evaldBool or endBool
            pass
        return endBool
        pass
    if operation.operation_type == "NOT":
        # Should always have only 1 child
        return not boolList[0]
        pass

def main(argv=None): # IGNORE:C0111
    '''Command line options.'''

    global proportion

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
        parser = argparse.ArgumentParser(prog='TGen-ASAP', description=program_license, formatter_class=argparse.RawTextHelpFormatter)
        required_group = parser.add_argument_group("required arguments")
        required_group.add_argument("-j", "--json", metavar="FILE", required=True, type=argparse.FileType('r'), help="JSON file of assay descriptions. [REQUIRED]")
        required_group.add_argument("-b", "--bam", metavar="FILE", required=True, type=argparse.FileType('rb'), default=sys.stdin, help="BAM file to analyze. [REQUIRED]")
        parser.add_argument("-d", "--depth", default=100, type=int, help="minimum read depth required to consider a position covered. [default: 100]")
        parser.add_argument("--breadth", default=0.8, type=float, help="minimum breadth of coverage required to consider an amplicon as present. [default: 0.8]")
        parser.add_argument("-p", "--proportion", type=float, help="minimum proportion required to call a mutation at a given locus. [default: 0.1]") #Don't explicitly set default because I need to be certain whether user set the value
        parser.add_argument("-m", "--mutation-depth", dest="mutdepth", default=5, type=int, help="minimum number of reads required to call a mutation at a given locus. [default: 5]")
        parser.add_argument("-V", "--version", action="version", version=program_version_message)
        parser.add_argument("-D", "--debug", action="store_true", default=False, help="write <sample_name>.log file with debugging information")
        parser.add_argument("-w", "--whole-genome", action="store_true", dest="wholegenome", default=False, help="JSON file uses a whole genome reference, so don't write out the consensus, depth, and proportion arrays for each sample")
        parser.add_argument("--allele-output-threshold", dest="allele_min_reads", default=8, type=int, help="cutoff of # of reads below which allels for amino acids and nucleotide alleles will not be output [default: 8]")
        parser.add_argument('-o', '--out', metavar="FILE", type=argparse.FileType('w'), default=sys.stdout, help="output filename [default: stdout]")
        parser.add_argument("--output-format", type=str.lower, choices=('xml', 'json'), default='xml', help="output format [default: xml]")
        parser.add_argument("--min-base-qual", dest="bqual", default=5, type=int, help="What is the minimum base quality score to use a position (phred scale, i.e. 10=90, 20=99, 30=99.9 accuracy) [default: 5]")
        parser.add_argument("--consensus-proportion", default=0.8, type=float, help="minimum proportion required to call at base at that position, else 'N'. [default: 0.8]")
        parser.add_argument("--fill-gaps", nargs="?", const="n", dest="gap_char", help="fill no coverage gaps in the consensus sequence [default: False], optional parameter is the character to use for filling [defaut: n]")
        parser.add_argument("--mark-deletions", nargs="?", const="_", dest="del_char", help="fill deletions in the consensus sequence [default: False], optional parameter is the character to use for filling [defaut: _]")

        # Process arguments
        args = parser.parse_args()

        json_fp = args.json
        bam_fp = args.bam
        bam_file = bam_fp.name
        depth = args.depth
        breadth = args.breadth
        proportion = args.proportion
        mutdepth = args.mutdepth
        debug = args.debug
        allele_min_reads = args.allele_min_reads
        wholegenome = args.wholegenome
        base_qual = args.bqual
        con_prop = args.consensus_proportion
        fill_gap_char = args.gap_char
        fill_del_char = args.del_char

        #out_dir = args.odir
        #if not out_dir:
        #    out_dir = os.getcwd()

        #out_dir = dispatcher.expandPath(out_dir)
        #if not os.path.exists(out_dir):
        #    os.makedirs(out_dir)

        operation_list = []
        operation_err = ""
        try:
            operation_list = assayInfo.parseOperation(args.json.name)
            pass
        except Exception as e:
            operation_err = str(e)
        assay_list = assayInfo.parseJSON(args.json.name)

        samdata = pysam.AlignmentFile(bam_fp.name, "rb")
        sample_dict = {}
        if 'RG' in samdata.header.to_dict() :
            sample_dict['name'] = samdata.header.to_dict()['RG'][0]['ID']
        else:
            sample_dict['name'] = os.path.splitext(os.path.basename(bam_fp.name))[0]
        sample_dict['mapped_reads'] = str(samdata.mapped)
        sample_dict['unmapped_reads'] = str(samdata.unmapped)
        sample_dict['unassigned_reads'] = str(samdata.nocoordinate)
        sample_dict['depth_filter'] = str(depth)
        sample_dict['proportion_filter'] = str(proportion)
        sample_dict['breadth_filter'] = str(breadth)
        sample_dict['mutation_depth_filter'] = str(mutdepth)
        # minidom.parseString will raise xml.parsers.expat.ExpatError: not well-formed (invalid token)
        # if json_file or bam_file contain the python string representation of a file-like object.
        # e.g. bam_file="<_io.BufferedReader name=\'/shared/Targeted_sequence_fastqs/ASAP/COD-10-24_S302.bam\'>"
        sample_dict['json_file'] = json_fp.name
        sample_dict['bam_file'] = bam_fp.name
        sample_node = ElementTree.Element("sample", sample_dict)

        if INFO or debug:
            if not os.path.isdir("./bamProcessorLogs"):
                os.mkdir("./bamProcessorLogs")
                pass
            logfile = "./bamProcessorLogs/%s.log" % sample_dict['name']
            pass

        if INFO:
            logging.basicConfig(level=logging.INFO,
                                format='%(asctime)s %(levelname)-8s %(message)s',
                                datefmt='%m/%d/%Y %H:%M:%S',
                                filename=logfile,
                                filemode='w')
            pass
        if debug:
            logging.basicConfig(level=logging.DEBUG,
                                format='%(asctime)s %(levelname)-8s %(message)s',
                                datefmt='%m/%d/%Y %H:%M:%S',
                                filename=logfile,
                                filemode='w')
        if operation_err != "":
            logging.info("Operation fetch did not succeed: "+operation_err)
            pass
        logging.info("----------------------bamProcessor STARTED---------------------------")
        logging.info("JSON: "+str(json_fp))
        logging.info("BAM: "+str(bam_fp.name))
        for assay in assay_list:
            assay_dict = {}
            assay_dict['name'] = assay.name
            assay_dict['type'] = assay.assay_type
            assay_dict['function'] = assay.target.function or ""
            assay_dict['gene'] = assay.target.gene_name or ""
            assay_dict['start'] = assay.target.start_position or ""
            assay_dict['end'] = assay.target.end_position or ""
            logging.info("Assay: "+str(assay.name))
            #offset is where the amplicon sits relative to a reference, subtract 1 to make 0-based position
            #report the lesser of start and end in case amplicon is on reverse strand
            try:
                offset = min(int(assay.target.start_position), int(assay.target.end_position))-1
                ref_positions = list(range(min(int(assay.target.start_position), int(assay.target.end_position)),max(int(assay.target.start_position), int(assay.target.end_position))+1))
                #remove the zero if present
                if 0 in ref_positions:
                    ref_positions.remove(0)
            except:
                offset = 0
                ref_positions = None
            assay_node = ElementTree.SubElement(sample_node, "assay", assay_dict)
            ref_name = assay.name
            reverse_comp = assay.target.reverse_comp
            for amplicon in assay.target.amplicons:
                logging.info("+++AMPLICON+++")
                if not ref_positions:
                    ref_positions = list(range(1, len(amplicon.sequence)+1))
                temp_file = None
                ref_name = assay.name + "_%s" % amplicon.variant_name if amplicon.variant_name else assay.name
                logging.info("Now Checking Amplicon: "+str(ref_name))
                amplicon_dict = {}
                seq_counter = None
                if samdata.closed:
                    samdata = pysam.AlignmentFile(bam_file, "rb")
                amplicon_dict['reads'] = str(samdata.count(ref_name))
                if amplicon.variant_name:
                    amplicon_dict['variant'] = amplicon.variant_name
                amplicon_node = ElementTree.SubElement(assay_node, "amplicon", amplicon_dict)
                if seq_counter:
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
                    # Warning: not designed to handle greater than 10 million X coverage
                    pileup = samdata.pileup(ref_name, max_depth=10000000, ignore_orphans=False, ignore_overlaps=False)
                    amplicon_data = _process_pileup(pileup, amplicon, depth, proportion, mutdepth, offset, wholegenome, base_qual, con_prop, fill_gap_char, fill_del_char)
                    if float(amplicon_data['breadth']) < breadth*100:
                        significance_node = amplicon_node.find("significance")
                        if significance_node is None:
                            significance_node = ElementTree.SubElement(amplicon_node, "significance")
                        if not significance_node.get("flag"):
                            significance_node.set("flag", "insufficient breadth of coverage")
                    # Handle SNPs
                    for snp in amplicon_data['SNPs']:
                        _add_snp_node(amplicon_node, snp)
                        # This would be helpful, but count_coverage is broken in python3 -- TODO: Revisit this
                        # print(samdata.count_coverage(ref_name, snp.position-1, snp.position))
                    del amplicon_data['SNPs']
                    _write_parameters(amplicon_node, amplicon_data)
                    if not wholegenome:
                        ref_positions_node = ElementTree.SubElement(amplicon_node, "ref_positions")
                        ref_positions_node.text = ",".join(str(n) for n in ref_positions)
                    # Handle ROIs
                    for roi in amplicon.ROIs:
                        roi_dict = _process_roi(roi, samdata, ref_name, smor, len(amplicon.sequence), reverse_comp)
                        _add_roi_node(amplicon_node, roi, roi_dict, depth, proportion, mutdepth, smor, offset, allele_min_reads)

                if temp_file and REMOVE_TEMP:
                    samdata.close()
                    os.remove(temp_file)
                    os.remove(temp_file+".bai")
                    samdata = pysam.AlignmentFile(bam_file, "rb")

        # Close File
        if samdata.is_open():
            samdata.close()

        # Handle Operations
        # For each true operation add a 'significance' into the output detailing the operation
        for operation in operation_list:
            if evaluateOperation(operation, sample_node):
                significance_node = ElementTree.SubElement(sample_node, "operation")
                significance_node.set("flag", str(operation))
                significance_node.set("message", str(operation.message))
                pass
            pass


        logging.info("Writing output: "+str(args.out))
        _write_output(args.out, sample_node, args.output_format)

    except KeyboardInterrupt:
        logging.info("Process ended via KeyboardInterrupt.")
        pass
    except Exception as e:
        if DEBUG or TESTRUN:
            raise(e)
        indent = len(cmdParser.program_name) * " "
        sys.stderr.write(cmdParser.program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        logging.info("An Exception Occured! "+str(e))
        return 2
    logging.info("-------------------------bamProcessor FINISHED EXIT(0)---------------------------")
    return 0

def _write_output(file_obj, xml_element, output_format='xml'):
    if output_format == 'xml':
        from xml.dom import minidom
        dom = minidom.parseString(ElementTree.tostring(xml_element))
        file_obj.write(dom.toprettyxml(indent="  "))
    elif output_format == 'json':
        xml_str = ElementTree.tostring(xml_element)
        # The 'sample' root node is discarded as an unnecessary layer for the JSON object.
        xml_obj = xmltodict.parse(xml_str)['sample']
        # FIXME: The output is en/decoded multiple times because it seemed
        # easier to use the json object_hook to ensure each key had a
        # a consistent type than to write a nested loop with type checks
        # and conversions modifying the object as it was traversed.
        #
        # Ideally the output should start as a python object that is
        # encoded to XML or JSON once.
        json_encoded_xml = json.loads(json.dumps(xml_obj), object_hook=cast_json_output_types)
        json.dump(json_encoded_xml, file_obj, separators=(',', ':'))
    else:
        raise Exception('unsupported output format: %s' % output_format)

# cast_json_output_types is a json decoder object_hook intended to be used on
# a ASAP output decoded from XML:
# - casts numbers from strings to float/int
# - keys that are expected to contain 0 to n elements are lists (or undefined)
#   eliminating the 1 element object case.
# - As a special addition, values that were stored in the XML as strings of
#   comma separated values are converted to an array of an appropriate type.
def cast_json_output_types(e):
    ## Sample
    if '@breadth_filter' in e:
        e['@breadth_filter'] = float(e['@breadth_filter'])
    if '@depth_filter' in e:
        e['@depth_filter'] = int(e['@depth_filter'])
    # @json_file
    if '@mapped_reads' in e:
        e['@mapped_reads'] = int(e['@mapped_reads'])
    # @name
    if '@proportion_filter' in e:
        e['@proportion_filter'] = float(e['@proportion_filter'])
    if '@unassigned_reads' in e:
        e['@unassigned_reads'] = int(e['@unassigned_reads'])
    if '@unmapped_reads' in e:
        e['@unmapped_reads'] = int(e['@unmapped_reads'])
    if 'assay' in e and not isinstance(e['assay'], list):
        e['assay'] = [(e['assay'])]

    ## Assay
    # @function
    # @gene
    # @name
    # @type
    # amplicon
    if 'amplicon' in e and not isinstance(e['amplicon'], list):
        e['amplicon'] = [e['amplicon']]

    ## Amplicon
    if '@reads' in e:
        e['@reads'] = int(e['@reads'])
    #if 'significance' in e and not isinstance(e['significance'], dict):
    # consensus_sequence
    if 'breadth' in e:
        e['breadth'] = float(e['breadth'])
    if 'depths' in e:
        e['depths'] = [int(v) for v in e['depths'].split(',')]
    if 'proportions' in e:
        e['proportions'] = [float(v) for v in e['proportions'].split(',')]
    if 'average_depth' in e:
        e['average_depth'] = float(e['average_depth'])
    if 'snp' in e and not isinstance(e['snp'], list):
        e['snp'] = [e['snp']]
    if 'region_of_interest' in e and not isinstance(e['region_of_interest'], list):
        e['region_of_interest'] = [e['region_of_interest']]

    ## SNP
    if '@depth' in e:
        e['@depth'] = int(e['@depth'])
    # @name
    if '@position' in e:
        e['@position'] = int(e['@position'])
    # @reference
    # snp_call
    if 'base_distribution' in e:
        e['base_distribution'] = {k: int(v) for k, v in e['base_distribution'].items()}

    ## SnpCall
    if '@count' in e:
        e['@count'] = int(e['@count'])
    if '@percent' in e:
        e['@percent'] = float(e['@percent'])
    # #text: "T"

    ## RegionOfInterest
    # TODO: what if {aa,nt}_sequence_distribution is set and None; should it default to an empty array? undefined? none?
    if e.get('aa_sequence_distribution'):
        e['aa_sequence_distribution'] = {k: int(v) for k, v in e['aa_sequence_distribution'].items()}
    if e.get('nt_sequence_distribution'):
        e['nt_sequence_distribution'] = {k: int(v) for k, v in e['nt_sequence_distribution'].items()}
    if '@changes' in e:
        e['@changes'] = int(e['@changes'])
    if 'mutation' in e and not isinstance(e['mutation'], list):
        e['mutation'] = [e['mutation']]

    ## Significance
    if 'significance' in e and isinstance(e['significance'], str):
        e['significance'] = {'#text': e['significance']}
    if '@resistance' in e:
        e['@resistance'] = e['@resistance'].split(',')

    return e

if __name__ == "__main__":
    if DEBUG:
        pass
    if TESTRUN:
        import doctest
        doctest.testmod()
    if PROFILE:
        import cProfile
        import pstats
        profile_filename = 'asap.newBamProcessor_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())

#!/usr/bin/env python3
# encoding: utf-8
'''
asap.bamProcessor -- Process BAM alignment files with an AssayInfo JSON file and generate XML for the results

asap.bamProcessor 

@author:     Darrin Lemmer

@copyright:  2015 TGen North. All rights reserved.

@license:    ACADEMIC AND RESEARCH LICENSE -- see ../LICENSE

@contact:    dlemmer@tgen.org
'''

import sys
import os
import re
import argparse
import logging

import pysam
from collections import Counter
from xml.etree import ElementTree
from skbio import DNA

import asap.dispatcher as dispatcher
import asap.assayInfo as assayInfo

__all__ = []
__version__ = 0.1
__date__ = '2015-07-16'
__updated__ = '2015-07-16'

DEBUG = 1
TESTRUN = 0
PROFILE = 0

def _write_parameters(node, data):
    for k, v in data.items():
        subnode = ElementTree.SubElement(node, k)
        subnode.text = str(v)
    return node

def _process_pileup(pileup, amplicon, depth, proportion):
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
            snp_call_proportion = alignment_call_proportion
        elif len(ordered_list) > 1:
            snp_call = ordered_list[1][0]
            snp_call_proportion = ordered_list[1][1] / pileupcolumn.n
        else:
            snp_call = snp_call_proportion = None
        consensus_seq += alignment_call if alignment_call_proportion >= proportion else "N"
        if position in snp_dict:
            for (name, reference, variant, significance) in snp_dict[position]:
                snp = {'name':name, 'position':str(position), 'depth':str(pileupcolumn.n), 'reference':reference, 'variant':variant, 'basecalls':base_counter}
                if base_counter[variant]/pileupcolumn.n >= proportion:
                    snp['significance'] = significance
                if not depth_passed:
                    snp['flag'] = "low coverage"
                snp_list.append(snp)
                #print("Found position of interest %d, reference: %s, distribution:%s" % (position, snp_dict[position][0], base_counter))
        elif depth_passed and snp_call and snp_call_proportion >= proportion:
            snp = {'name':'unknown', 'position':str(position), 'depth':str(pileupcolumn.n), 'reference':reference_call, 'variant':snp_call, 'basecalls':base_counter}
            if 0 in snp_dict:
                (name, *rest, significance) = snp_dict[0][0]
                snp['name'] = name
                snp['significance'] = significance
            snp_list.append(snp)
            #print("SNP found at position %d: %s->%s" % (position, reference_call, alignment_call))
    
    pileup_dict['consensus_sequence'] = consensus_seq
    pileup_dict['breadth'] = str(breadth_positions/amplicon_length * 100)
    pileup_dict['depths'] = ",".join(depth_array)
    pileup_dict['proportions'] = ",".join(prop_array)
    pileup_dict['SNPs'] = snp_list
    pileup_dict['average_depth'] = str(avg_depth_total/avg_depth_positions) if avg_depth_positions else "0"
    return pileup_dict 

def _write_xml(root, xml_file):
    from xml.dom import minidom
    #print(ElementTree.dump(root))
    dom = minidom.parseString(ElementTree.tostring(root))
    output = open(xml_file, 'w')
    output.write(dom.toprettyxml(indent="    "))
    output.close()
    return xml_file

def _create_snp_dict(amplicon):
    snp_dict = {}
    for snp in amplicon.SNPs:
        name = snp.name if snp.name else "position of interest"
        if snp.position in snp_dict:
            snp_dict[snp.position].append((name, snp.reference, snp.variant, snp.significance))
        else:
            snp_dict[snp.position] = [(name, snp.reference, snp.variant, snp.significance)]
    return snp_dict

def _add_snp_node(parent, snp):
    snp_attributes = {k:snp[k] for k in ('name', 'position', 'depth', 'reference')}
    snp_node = ElementTree.SubElement(parent, 'snp', snp_attributes)
    base_counter = snp['basecalls']
    snpcall = snp['variant']
    depth = int(snp['depth'])
    snpcount = base_counter[snpcall]
    snpcall_node = ElementTree.SubElement(snp_node, 'snp_call', {'count':str(snpcount), 'percent':str(snpcount/depth*100)})
    snpcall_node.text = snpcall
    if 'significance' in snp or 'flag' in snp:
        significance_node = ElementTree.SubElement(snp_node, 'significance')
        if 'significance' in snp:
            significance_node.text = snp['significance'].message
            if snp['significance'].resistance:
                significance_node.set("resistance", snp['significance'].resistance)
        if 'flag' in snp:
            significance_node.set('flag', snp['flag'])
    ElementTree.SubElement(snp_node, 'base_distribution', {k:str(v) for k,v in base_counter.items()})
    return snp_node

def _process_roi(roi, samdata, amplicon_ref, reverse_comp=False):
    roi_dict = {'region':roi.position_range}
    range_match = re.search('(\d*)-(\d*)', roi.position_range)
    if not range_match:
        return roi_dict
    start = int(range_match.group(1)) - 1
    end = int(range_match.group(2))
    expected_length = end - start
    aa_sequence_counter = Counter()
    nt_sequence_counter = Counter()
    depth = 0
    for read in samdata.fetch(amplicon_ref, start, end):
        rstart = read.reference_start
        alignment_length = read.get_overlap(start, end)
        #throw out reads that either have gaps in the ROI or don't cover the whole ROI
        if alignment_length != expected_length:
            continue
        if rstart <= start:
            for (qpos, rpos) in read.get_aligned_pairs():
                if rpos == start:
                    qstart = qpos
                if rpos == end:
                    qend = qpos
            #throw out reads with insertions in the ROI
            if not qend or not qstart or qend-qstart != expected_length:
                continue
            nt_sequence = DNA(read.query_alignment_sequence[qstart:qend])
            if reverse_comp:
                nt_sequence = nt_sequence.reverse_complement()
            #scikit-bio doesn't support translating degenerate bases currently, so we will just throw out reads with degenerates for now
            if nt_sequence.has_degenerates(): 
                continue
            aa_sequence = nt_sequence.translate()
            aa_string = str(aa_sequence).replace('*', 'x')
            if aa_string:
                nt_sequence_counter.update([str(nt_sequence)])
                aa_sequence_counter.update([aa_string])
                depth += 1
    if len(aa_sequence_counter) == 0:
        roi_dict['flag'] = "region not found"
        return roi_dict
    aa_consensus = aa_sequence_counter.most_common(1)[0][0]
    nt_consensus = nt_sequence_counter.most_common(1)[0][0]
    num_changes = 0
    reference = roi.aa_sequence
    consensus = aa_consensus
    if roi.nt_sequence:
        reference = roi.nt_sequence
        consensus = nt_consensus
    for i in range(len(reference)):
        if len(consensus) <= i or reference[i] != consensus[i]:
            num_changes += 1
    roi_dict['most_common_aa_sequence'] = aa_consensus
    roi_dict['most_common_nt_sequence'] = nt_consensus
    roi_dict['reference'] = reference
    roi_dict['changes'] = str(num_changes)
    roi_dict['aa_sequence_distribution'] = aa_sequence_counter
    roi_dict['nt_sequence_distribution'] = nt_sequence_counter
    roi_dict['depth'] = str(depth)
    return roi_dict

def _add_roi_node(parent, roi, roi_dict, depth, proportion):
    if "flag" in roi_dict:
        roi_node = ElementTree.SubElement(parent, "region_of_interest", {'region':roi_dict['region']})
        significance_node = ElementTree.SubElement(roi_node, "significance", {'flag':roi_dict['flag']})
        if roi.significance.resistance:
            significance_node.set("resistance", roi.significance.resistance)
        return roi_node
    roi_attributes = {k:roi_dict[k] for k in ('region', 'reference', 'depth')}
    roi_node = ElementTree.SubElement(parent, "region_of_interest", roi_attributes)
    aa_seq_counter = roi_dict['aa_sequence_distribution']
    aa_seq_count = aa_seq_counter[roi_dict['most_common_aa_sequence']]
    aa_seq_node = ElementTree.SubElement(roi_node, "amino_acid_sequence", {'count':str(aa_seq_count), 'percent':str(aa_seq_count/int(roi_dict['depth'])*100)})
    aa_seq_node.text = roi_dict['most_common_aa_sequence']
    nt_seq_counter = roi_dict['nt_sequence_distribution']
    nt_seq_count = nt_seq_counter[roi_dict['most_common_nt_sequence']]
    nt_seq_node = ElementTree.SubElement(roi_node, "nucleotide_sequence", {'count':str(nt_seq_count), 'percent':str(nt_seq_count/int(roi_dict['depth'])*100)})
    nt_seq_node.text = roi_dict['most_common_nt_sequence']
    significant = False
    for mutation in roi.mutations:
        if roi.nt_sequence:
            count = nt_seq_counter[mutation]
            mutant_proportion = count/int(roi_dict['depth'])
        else:
            count = aa_seq_counter[mutation]
            mutant_proportion = count/int(roi_dict['depth'])
        mutation_node = ElementTree.SubElement(roi_node, 'mutation', {'name':str(roi.name)+mutation, 'count':str(count), 'percent':str(mutant_proportion*100)})
        mutation_node.text = mutation
        if mutant_proportion >= proportion:
            significant = True
    if significant:
        significance_node = ElementTree.SubElement(roi_node, "significance")
        significance_node.text = roi.significance.message
        if roi.significance.resistance:
            significance_node.set("resistance", roi.significance.resistance)
        if int(roi_dict['depth']) < depth:
            significance_node.set("flag", "low coverage")
    elif 'changes' in roi_dict and int(roi_dict['changes']) > 0:
        significance_node = ElementTree.SubElement(roi_node, "significance", {'changes':roi_dict['changes']})
        significance_node.text = roi.significance.message                       
        if roi.significance.resistance:
            significance_node.set("resistance", roi.significance.resistance)
    ElementTree.SubElement(roi_node, 'aa_sequence_distribution', {k:str(v) for k,v in aa_seq_counter.items()})
    ElementTree.SubElement(roi_node, 'nt_sequence_distribution', {k:str(v) for k,v in nt_seq_counter.items()})
    return roi_node

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
        parser.add_argument("-p", "--proportion", default=0.1, type=float, help="minimum proportion required to call a SNP at a given position. [default: 0.1]")
        parser.add_argument("-V", "--version", action="version", version=program_version_message)
     
        # Process arguments
        args = parser.parse_args()

        json_fp = args.json
        bam_fp = args.bam
        out_fp = args.out
        depth = args.depth
        breadth = args.breadth
        proportion = args.proportion
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
        if 'RG' in samdata.header :
            sample_dict['name'] = samdata.header['RG'][0]['ID']
        else:
            sample_dict['name'] = os.path.splitext(os.path.basename(bam_fp))[0]
        sample_dict['mapped_reads'] = str(samdata.mapped)
        sample_dict['unmapped_reads'] = str(samdata.unmapped)
        sample_dict['unassigned_reads'] = str(samdata.nocoordinate)
        sample_node = ElementTree.Element("sample", sample_dict)

        #out_fp = os.path.join(out_dir, sample_dict['name']+".xml")
        
        for assay in assay_list:
            assay_dict = {}
            assay_dict['name'] = assay.name
            assay_dict['type'] = assay.assay_type
            assay_node = ElementTree.SubElement(sample_node, "assay", assay_dict)
            ref_name = assay.name
            reverse_comp = assay.target.reverse_comp
            for amplicon in assay.target.amplicons:
                ref_name = assay.name + "_%s" % amplicon.variant_name if amplicon.variant_name else assay.name
                amplicon_dict = {}
                amplicon_dict['reads'] = str(samdata.count(ref_name))
                if amplicon.variant_name:
                    amplicon_dict['variant'] = amplicon.variant_name
                amplicon_node = ElementTree.SubElement(assay_node, "amplicon", amplicon_dict)
                if samdata.count(ref_name) == 0:
                    significance_node = ElementTree.SubElement(amplicon_node, "significance", {"flag":"no coverage"})
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
                    amplicon_data = _process_pileup(pileup, amplicon, depth, proportion)
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
                        roi_dict = _process_roi(roi, samdata, ref_name, reverse_comp)
                        _add_roi_node(amplicon_node, roi, roi_dict, depth, proportion)

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

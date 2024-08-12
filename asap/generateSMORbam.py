#!/usr/bin/env python3
# encoding: utf-8
'''
asap.generateSMORbam -- Generate a new bam file for a SMOR analysis, containing the SMOR consensus reads

asap.generateSMORbam 

@author:     Darrin Lemmer

@copyright:  2024 TGen North. All rights reserved.

@license:    ACADEMIC AND RESEARCH LICENSE -- see ../LICENSE

@contact:    dlemmer@tgen.org
'''

import sys
import os
import re
import argparse
import logging
import pysam
from operator import attrgetter
from statistics import mode

from asap import __version__ 

__all__ = []
__date__ = '2024-07-11'
__updated__ = '2024-07-11'

DEBUG = 1
TESTRUN = 0
PROFILE = 0

def pairwise(iterable):
    from itertools import tee
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

def _find_overlap_region(reads):
    start_list = []
    end_list = []
    for read in reads:
        start_list.append(read.reference_start)
        end_list.append(read.reference_end)
    if start_list and end_list:
        return (mode(start_list), mode(end_list))
    else:
        return (None, None)

#Not used
def _get_consensus2(seq1, seq2, qual1, qual2, fill_char, base_qual):
    consensus = ""
    quals = []

    for i in range(min(len(seq1), len(seq2))):
        b1 = seq1[i]
        b2 = seq2[i]
        q1 = qual1[i]
        q2 = qual2[i]
        if b1 == b2 and q1 >= base_qual and q2 >= base_qual:
            consensus = consensus + b1
        else:
            consensus = consensus + fill_char
        quals.append(min(q1, q2))
                    
    return (consensus, quals)

def _get_consensus(read, pair, start, end, fill_char):
    consensus = ""
    quals = []
    cigar = []
    cigar_length = 0
    cigar_op = 0
    end = min(end, read.reference_end, pair.reference_end)
    read1_alignment = read.get_aligned_pairs()
    read2_alignment = pair.get_aligned_pairs()
    read1_iter = iter(read1_alignment)
    read2_iter = iter(read2_alignment)
    read1_ptr = next(read1_iter)
    read2_ptr = next(read2_iter)
    while read1_ptr[1] == None or read1_ptr[1] < start or read1_ptr[1] < read2_ptr[1]:
        read1_ptr = next(read1_iter, None)
    while read2_ptr[1] == None or read2_ptr[1] < start or read2_ptr[1] < read1_ptr[1]:
        read2_ptr = next(read2_iter, None)

    if read1_ptr[1] > start or read2_ptr[1] > start: #We are starting with a deletion
        cigar_op = 2
        cigar_length += (max(read1_ptr[1], read2_ptr[1]) - start)

    #Both pointers /should/ be at the overlap start now
    while read1_ptr[1] != read2_ptr[1]: # deletion in one of the reads
        consensus += fill_char
        quals.append(0)
        if cigar_op != 0: #changing op, update cigar tuples and start a new one
            cigar.append((cigar_op, cigar_length))
            cigar_op = 0
            cigar_length = 0
        cigar_length += 1
        if read1_ptr[0] == None:
            read1_ptr = next(read1_iter)
        if read2_ptr[0] == None:
            read2_ptr = next(read2_iter)

    #Both pointers /should/ be at the same reference position now
    while read1_ptr[1] == read2_ptr[1]:
        if read1_ptr[0] == None and read2_ptr[0] == None: #Deletion in both, skip
            if cigar_op != 2: #changing op, update cigar tuples and start a new one
                cigar.append((cigar_op, cigar_length))
                cigar_op = 2
                cigar_length = 0
        elif read1_ptr[0] == None and read2_ptr[0] != None: #Deletion only in read1
            consensus += fill_char
            quals.append(pair.query_qualities[read2_ptr[0]])
            if cigar_op != 0: #changing op, update cigar tuples and start a new one
                cigar.append((cigar_op, cigar_length))
                cigar_op = 0
                cigar_length = 0
        elif read1_ptr[0] != None and read2_ptr[0] == None: #Deletion only in read2
            consensus += fill_char
            quals.append(read.query_qualities[read1_ptr[0]])
            if cigar_op != 0: #changing op, update cigar tuples and start a new one
                cigar.append((cigar_op, cigar_length))
                cigar_op = 0
                cigar_length = 0
        else: #Check for matching sequence
            if read1_ptr[1] == None: #Insertion in both
                if cigar_op != 1: #changing op, update cigar tuples and start a new one
                    cigar.append((cigar_op, cigar_length))
                    cigar_op = 1
                    cigar_length = 0
            else: #Match in both
                if cigar_op != 0: #changing op, update cigar tuples and start a new one
                    cigar.append((cigar_op, cigar_length))
                    cigar_op = 0
                    cigar_length = 0
            read1_base = read.query_sequence[read1_ptr[0]]
            read2_base = pair.query_sequence[read2_ptr[0]]
            if read1_base == read2_base:
                consensus += read1_base
            else:
                consensus += fill_char
            quals.append(min(read.query_qualities[read1_ptr[0]],pair.query_qualities[read2_ptr[0]]))
        cigar_length += 1
        #Advance the pointers to the next position
        if read1_ptr[1] == None or (read1_ptr[1] < end-1):# and read1_ptr[0] != None and read1_ptr[0] < read.query_alignment_end):
            read1_ptr = next(read1_iter, None)
        if read2_ptr[1] == None or read2_ptr[1] < end-1:
            read2_ptr = next(read2_iter, None)
        while read1_ptr[1] != read2_ptr[1]: # insertion in one of the reads, skip ahead until we are in sync
            if read1_ptr[1] == None:
                read1_ptr = next(read1_iter, None)
            if read2_ptr[1] == None:
                read2_ptr = next(read2_iter, None)
        if read1_ptr == None or read2_ptr == None or (read1_ptr[1] == end-1 and read2_ptr[1] == end-1):
            break

    cigar.append((cigar_op, cigar_length))
    return (consensus, quals, cigar)

#Not used
def _generate_cigar(offset, length, orig_cigar):
    new_cigar = []
    temp_cigar = []
    for element in orig_cigar:
        element_length = element[1]
        # Check if the length is less than or equal to the remaining offset
        if element_length <= offset:
            # Subtract the length from the offset and remove the tuple element
            if element[0] != 2:
                offset -= element_length
        else:
            # If the length is larger, subtract the remaining offset and update the element
            new_element = (element[0], element_length - offset)
            temp_cigar.append(new_element)
            # Add the rest of the list to the result and break the loop
            temp_cigar.extend(orig_cigar[orig_cigar.index(element) + 1:])
            break
    for element in temp_cigar:
        element_length = element[1]
        # Check if the length is less than or equal to the remaining offset
        if element_length <= length:
            # add in elements until we've covered the alignment length
            new_cigar.append(element)
            if element[0] != 2:
                length -= element_length
        else:
            # If the length is larger, subtract the difference and update the element
            new_element = (element[0], element_length - (element_length - length))
            new_cigar.append(new_element)
            # Break the loop to ignore the rest of the original cigar
            break
    return (new_cigar)   

def _write_bam(samdata, out_file, fill_char, base_qual):
    outdata = pysam.AlignmentFile(out_file, "wb", template=samdata)
    
    for ref_name in samdata.references:
        (start, end) = _find_overlap_region(samdata.fetch(ref_name))
        reads = iter(sorted(samdata.fetch(ref_name), key=attrgetter('query_name')))
        for read, pair in pairwise(reads):
            if read.query_name != pair.query_name:
                continue
            if read.reference_end == None or pair.reference_end == None:
                continue
            if read.reference_end < pair.reference_start or read.reference_start > pair.reference_end:
                continue
            read1_alignment = read.get_aligned_pairs(True)
            read2_alignment = pair.get_aligned_pairs(True)
            if not read1_alignment or not read2_alignment:
                continue
            read1_start = 0
            read1_end = 0
            read2_start = 0
            read2_end = 0
            for match in read1_alignment:
                if match[1] == start:
                    read1_start = match[0]
                if match[1] == end-1:
                    read1_end = match[0]
            for match in read2_alignment:
                if match[1] == start:
                    read2_start = match[0]
                if match[1] == end-1:
                    read2_end = match[0]
            read_seq = read.query_sequence[read1_start:read1_end]
            pair_seq = pair.query_sequence[read2_start:read2_end]
            read_qual = read.query_qualities[read1_start:read1_end]
            pair_qual = pair.query_qualities[read2_start:read2_end]

            (consensus, quals, cigar) = _get_consensus(read, pair, start, end, fill_char)
            #(consensus, quals) = _get_consensus2(read_seq, pair_seq, read_qual, pair_qual, fill_char, base_qual)
            #if read.query_length + read.get_cigar_stats()[0][2] <= pair.query_length + pair.get_cigar_stats()[0][2]:
            #    new_cigar = _generate_cigar(read1_start, len(consensus), read.cigartuples)
            #else:
            #    new_cigar = _generate_cigar(read2_start, len(consensus), pair.cigartuples)
            new_read = pysam.AlignedSegment()
            new_read.is_paired = False
            new_read.query_name = read.query_name
            new_read.reference_id = read.reference_id
            new_read.reference_start = start
            new_read.query_sequence = consensus
            new_read.cigartuples = cigar
            new_read.query_qualities = quals
            outdata.write(new_read)
    outdata.close()
    pysam.sort("-o", out_file, out_file)
    pysam.index(out_file)
    return (out_file)

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
        parser = argparse.ArgumentParser(description=program_license, formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument("-b", "--bam", metavar="FILE", required=True, type=argparse.FileType('rb'), help="bam file to process. [REQUIRED]")
        parser.add_argument("-o", "--out", metavar="FILE", help="new bam file to write. [default: ./{orig_bam}_SMOR.bam]")
        parser.add_argument("-c", "--fill-character", default="N", type=str, dest="fill", help="character to use for overlap positions that don't match [default: N]")
        parser.add_argument("-q", "--min-base-qual", dest="bqual", default=0, type=int, help="minimum base quality score to use a position in each read (Phred scale, i.e. 10=90, 20=99, 30=99.9 percent accuracy) [default: 0]")
        
        parser.add_argument('-V', '--version', action='version', version=program_version_message)

        # Process arguments
        args = parser.parse_args()

        bam_fp = args.bam
        out_file = args.out
        fill_char = args.fill
        base_qual = args.bqual

        samdata = pysam.AlignmentFile(bam_fp.name, "rb")

        if not out_file:
            out_file = "%s_SMOR.bam" % (os.path.splitext(os.path.basename(samdata.filename.decode("utf-8")))[0])
            
        _write_bam(samdata, out_file, fill_char, base_qual)

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
        profile_filename = 'asap.generateSMORbam_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    main()

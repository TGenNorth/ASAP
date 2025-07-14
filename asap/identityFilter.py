#!/usr/bin/env python3
# encoding: utf-8
'''
asap.identityFilter -- Generate a new bam file just like the old one, but with primer sequences masked

asap.identityFilter 

@author:     Darrin Lemmer

@copyright:  2025 TGen North. All rights reserved.

@license:    ACADEMIC AND RESEARCH LICENSE -- see ../LICENSE

@contact:    dlemmer@tgen.org
'''

import sys
import os
import re
import argparse
import logging
import pysam
import numpy as np
import array as arr
from operator import attrgetter
from statistics import mode

from asap import __version__ 

__all__ = []
__date__ = '2025-03-21'
__updated__ = '2025-03-21'

DEBUG = 1
TESTRUN = 0
PROFILE = 0

def pairwise(iterable):
    from itertools import tee
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

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
    #This needs more work. Is it worth it?

def _mark_read_unaligned(read):
    read.is_unmapped = True
    read.is_secondary = False
    read.is_supplementary = False
    read.reference_id = -1
    read.reference_start = -1
    read.cigar = []
    read.next_reference_id = -1
    read.next_reference_start = -1
    read.mapping_quality = 0
    read.tags = [tag for tag in read.tags if tag[0] != 'NM' and tag[0] != 'MD' and tag[0] != 'AS']
    return(read)

def _identity_filter(samdata, ref_names, percid, merge, out_fp):
    outdata = pysam.AlignmentFile(out_fp, "wb", template=samdata)
    discarded_reads = 0
    #seq_counter = Counter()
    aligned_reads = []
    
    # if user didn't specify any specific refs, apply percid filter to all
    if ref_names == None:
        ref_names = samdata.references
    
    # Iterate through all the references that need to be filtered
    for read in samdata.fetch(until_eof=True):
        if read.is_unmapped:
            logging.info("Read %s is unmapped -- copying it over, as is...." % read.query_name);
            outdata.write(read)
            continue
        if read.reference_name in ref_names:    
            length = read.infer_query_length(False)
            logging.info("Checking %s against reference %s" % (read.query_name, read.reference_name))
            logging.info("\tAligned length %i, total read length %i" % (read.query_alignment_length or -1, length or -1))
            if read.query_alignment_length / length >= percid: #Quick check that the aligned length even passes threshold
                matches = 0
                for (qpos, rpos, seq) in read.get_aligned_pairs(with_seq=True):
                    query = read.query_sequence[qpos] if qpos else "None"
                    #logging.info("\tqpos: %i\trpos: %i\tseq: %s\tquery[qpos]: %s" % (qpos or -1, rpos or -1, seq, query))
                    #if there is a gap in the alignment, extend the length of the query or reference accordingly
                    if rpos is None:
                        pass #amp_length += 1
                    elif qpos is None:
                        length += 1
                    else:
                        if read.query_sequence[qpos].upper() == seq.upper():
                            matches += 1
                if matches / length >= percid: #Using length instead of amp_length to compare to query instead of reference
                    logging.info("\t\tFound %i matches out of %i, keeping..." % (matches, length))
                    outdata.write(read)
                else:
                    logging.info("\t\tFound %i matches out of %i, marking as unaligned..." % (matches, length))
                    discarded_reads += 1
                    #seq_counter.update([read.query_sequence])
                    outdata.write(_mark_read_unaligned(read))
            else: #aligned proportion below threshold
                logging.info("\t\tAlignment too short, marking as unaligned...")
                discarded_reads += 1
                #seq_counter.update([read.query_sequence])
                outdata.write(_mark_read_unaligned(read))
        else: # Read is not aligned to a reference we are verifying, let it go
            logging.info("Read %s is aligned to a reference we aren't checking -- copying it over, as is...." % read.query_name);
            outdata.write(read)
    outdata.close()
    return (outdata, discarded_reads)

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
        parser.add_argument("-i", "--identity", metavar="float", dest = "percid", required=True, help="minimum percent identity required to keep aligned read. [REQUIRED]")
        parser.add_argument("-r", "--ref", dest="ref_names", action="append", nargs="+", help="name of the reference contig(s) for which identity is calculated; if omitted, apply to all contigs. May be specified multiple times")
        #parser.add_argument("-m", "--merge", action="store_true", default=False, help="merge paired reads before calculating identity. [default: False]")
        parser.add_argument("-o", "--out", metavar="FILE", help="new bam file to write. [default: ./{orig_bam}_identityFiltered.bam]")
        parser.add_argument('-V', '--version', action='version', version=program_version_message)

        # Process arguments
        args = parser.parse_args()

        bam_fp = args.bam
        out_fp = args.out
        percid = float(args.percid)
        ref_names = args.ref_names
        merge = False #args.merge

        logfile = "identity_filtering.log"
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s %(levelname)-8s %(message)s',
                            datefmt='%m/%d/%Y %H:%M:%S',
                            filename=logfile,
                            filemode='w')

        samdata = pysam.AlignmentFile(bam_fp.name, "rb")

        if not out_fp:
            out_fp = "%s_identityFiltered.bam" % (os.path.splitext(os.path.basename(samdata.filename.decode("utf-8")))[0])
     
        (samout, discarded_reads) = _identity_filter(samdata, ref_names, percid, merge, out_fp)
                
        #bam_file_out_sorted = "%s_sorted.bam" % (os.path.splitext(os.path.basename(samout.filename.decode("utf-8")))[0])
        
        pysam.sort("-o", out_fp, out_fp)
        pysam.index(out_fp)

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
        profile_filename = 'asap.identityFilter_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    main()

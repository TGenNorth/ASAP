#!/usr/bin/env python3
# encoding: utf-8
'''
asap.maskPrimers -- Generate a new bam file just like the old one, but with primer sequences masked

asap.maskPrimers 

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

def _primer_mask(samdata, primer_file, wiggle, mask_bases, ponlybam, outfile):
    # TODO: Ideally for smor it should use the same primer pair, so will have to match reads and then search, will at this at some stage
    # assumptions
    # you know that primerF is on read 1 and primerR is for read 2 (as you added the adapters like this)
    # you want to keep singletons (these could be easily removed later)
    logging.info("Starting primer_mask function")
    outdata = pysam.AlignmentFile(outfile, "wb", template=samdata)
    out = open('primer_masking.tsv', 'w')
    out.write('refname\tread_id\tprimer_name\tprimer_region\tmasked_sequence\n')
    
    # process primer file - add some error handling here
    try:
        primers = np.loadtxt(str(primer_file), delimiter="\t",
          dtype={'names': ('CHROM', 'PrimerName', 'PrimerDirection', 'Start', 'End'),
          'formats': ('<U100', '<U100', 'U1', 'int', 'int')}, skiprows=1)
    except ValueError:
        logging.error("Incorrect primer file format")
        return samdata
    primers["PrimerDirection"] = np.char.upper(primers["PrimerDirection"])
    primer_stats = []
    # TODO add ponlybam option to only emit reads with primer sequence, will need to deal with pairs in this case
    # for each ref in bam
    for chrom in samdata.references:
        # check that all chroms are accounted for in input file
        if chrom in primers["CHROM"]:
            # For each primer set
            #from inputus--------
            # primers["CHROM"][primers["CHROM"] == chrom]
            forward_primer_set_strt = primers["Start"][np.multiply(primers["CHROM"] == chrom, primers["PrimerDirection"] == "F")]
            forward_primer_set_end = primers["End"][np.multiply(primers["CHROM"] == chrom, primers["PrimerDirection"] == "F")]
            reverse_primer_set_strt = primers["Start"][np.multiply(primers["CHROM"] == chrom, primers["PrimerDirection"] == "R")]
            reverse_primer_set_end = primers["End"][np.multiply(primers["CHROM"] == chrom, primers["PrimerDirection"] == "R")]
            #from inputus end--------
            forward_primer_set_strt = forward_primer_set_strt - wiggle
            reverse_primer_set_end = reverse_primer_set_end + wiggle
            # make sure no negative
            forward_primer_set_strt[forward_primer_set_strt < 0] = 0
            #make sure not longer than reference
            ref_len = samdata.get_reference_length(chrom)
            forward_primer_set_end[forward_primer_set_end > ref_len] = ref_len
            reverse_primer_set_end[reverse_primer_set_end > ref_len] = ref_len
            no_primer = 0
            primer_found = 0
            for read in samdata.fetch(chrom, until_eof=True):
                try: # It seems this happens when only one read in the pair is aligned, then the mate is still associated with the CHROM but has no alignemnt position, causing min() to fail as no value in it
                    # read.query_alignment_start (ead.query_alignment_end) is what base of the read is the first thats aligned to the reference
                    align_start = min(read.get_reference_positions()) #+ read.query_alignment_start #first base of read that is aligned, might be useful if we consider that adapters have been remove, not using now
                    align_end = max(read.get_reference_positions())
                except Exception as e:
                    no_primer += 1
                    # read.query_qualities = arr.array("B", [0] * len(read.query_qualities))
                    # read.query_sequence = "N" * len(read.query_sequence)
                    outdata.write(read) # this is here as to keep pairs, the above should probably be added?
                    out.write(f'{chrom}\t{read.query_name}\tNone\t\t{read.query_sequence}\n')
                    continue
                if read.is_read1:
                    # If read aligns within 'wiggle' nts of primer sequence
                    read_start_in_primer = np.multiply(align_start >= forward_primer_set_strt, align_start <= forward_primer_set_end)
                    if any(read_start_in_primer):
                        primer_found += 1
                        # deal with if more than one is true using the max
                        # to deal with if whole read in in the primer section, this can happen if short reads not filtered or if incorrect primer file, and not removing as then need to remove mate; TODO should output stats to log
                        # (align_end if int(forward_primer_set_end[tmp_boo].max()) > align_end else int(forward_primer_set_end[tmp_boo].max())) - align_start + read.query_alignment_start
                        # need to get the rightmost query position aligned to the primer area
                        primer_end_ref_pos = int(forward_primer_set_end[read_start_in_primer].max())

                        # Find the index where reference position equals the end of the primer, to get the corresponding query position
                        # if there is in indel there (read position = None), we may need to adjust the position
                        aligned_pairs = read.get_aligned_pairs()
                        target_idx = next((i for i, align in enumerate(aligned_pairs) if align[1] == primer_end_ref_pos), None)

                        if target_idx is not None:
                            # Work backwards from that index to find first non-None query position
                            mask_end = next(
                                (aligned_pairs[i][0] for i in range(target_idx, -1, -1) if aligned_pairs[i][0] is not None), align_end
                            )
                        else: # read doesn't align all the way to end of primer, so just go to the end of the read alignment
                            mask_end = align_end

                        mask_end = len(read.query_sequence) if mask_end > len(read.query_sequence) else mask_end
                        # This will work if using qual later for calling
                        read.query_qualities[:mask_end] = arr.array("B", [0] * mask_end)
                        if mask_bases:
                            qual_store = read.query_qualities
                            read.query_sequence = "N" * len(read.query_sequence[:mask_end]) + read.query_sequence[mask_end:]
                            read.query_qualities = qual_store
                        out.write(f'{chrom}\t{read.query_name}\tPrimerName\t0:{mask_end}\t{read.query_sequence}\n')
                    else:
                        # TODO
                        # else, mark as a fail. If fails > X% of reads then retry with larger wiggle?
                        no_primer += 1
                        out.write(f'{chrom}\t{read.query_name}\tNone\t\t{read.query_sequence}\n')
                elif read.is_read2:
                    read_end_in_primer = np.multiply(align_end >= reverse_primer_set_strt, align_end <= reverse_primer_set_end)
                    if any(read_end_in_primer):
                        primer_found += 1
                        # deal with if more than one is true using the min
                        primer_start_ref_pos = int(reverse_primer_set_strt[read_end_in_primer].min())
                        
                        # Find the index where reference position equals the start of the primer, to get the corresponding query position
                        # if there is in indel there (read position = None), we may need to adjust the position
                        aligned_pairs = read.get_aligned_pairs()
                        target_idx = next((i for i, align in enumerate(aligned_pairs) if align[1] == primer_start_ref_pos), None)

                        if target_idx is not None:
                            # Work forwards from that index to find first non-None query position
                            mask_start = next(
                                (aligned_pairs[i][0] for i in range(target_idx, len(aligned_pairs)) if aligned_pairs[i][0] is not None), align_start
                            )
                        else: # read doesn't align all the way to start of primer, so just go to the start of the read alignment
                            mask_start = align_start

                        read.query_qualities[mask_start:read.query_length] = arr.array("B", [0] * len(read.query_qualities[mask_start:read.query_length]))
                        if mask_bases:
                            qual_store = read.query_qualities
                            read.query_sequence = read.query_sequence[:mask_start] + "N" * len(read.query_sequence[mask_start:])
                            try:
                                read.query_qualities = qual_store
                            except ValueError as e:
                        out.write(f'{chrom}\t{read.query_name}\tPrimerName\t{mask_start}:{read.query_length}\t{read.query_sequence}\n')
                    else:
                        no_primer += 1
                        out.write(f'{chrom}\t{read.query_name}\tNone\t\t{read.query_sequence}\n')
                else:
                    logging.debug("Aberrant read: %s" % read.query_name)
                outdata.write(read)
            primer_stats.append([chrom, primer_found, no_primer])
        else:
            logging.info("No primers found for: %s" % chrom)
            for read in samdata.fetch(chrom, until_eof=True):
                outdata.write(read)
    logging.info("CHROM, Primer Found, Primer Missing")
    logging.info(primer_stats)
    outdata.close() #only aligned reads
    samdata.close()
    if mask_bases: # need to sort as trimming may have changed coordinates
        bam_file_out_sorted = "%s_primerMasked_sorted.bam" % (os.path.splitext(os.path.basename(samdata.filename.decode("utf-8")))[0])
        pysam.sort("-o", bam_file_out_sorted, outfile)
        outfile = bam_file_out_sorted
    pysam.index(outfile)
    logging.info("Wrote primer masking alignment only file: %s" % outfile)
    return outfile

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
        parser.add_argument("-p", "--primer-file", metavar="FILE", dest = "primers", help="primer file to use for primer masking. [REQUIRED]")
        parser.add_argument("-o", "--out", metavar="FILE", help="new bam file to write. [default: ./{orig_bam}_primerMasked.bam]")
        parser.add_argument("--wiggle", dest="wiggle", default=9, type=int, help="How many nucleotides outside the primer window should be used to identify primer sequences [default: 9]")
        parser.add_argument("--mask-bam", dest="maskbam", action="store_true", default=True, help="change primer sequences in the alignment file to 'Ns' [default]")
        parser.add_argument("--no-mask-bam", dest="maskbam", action="store_false", help="don't modify primer sequences in the alignment file")
        parser.add_argument("--primer-only", dest="primeronly", action="store_true", default=False, help="only keep sequences with primers")
        parser.add_argument("--no-primer-only", dest="primeronly", action="store_false", default=True, help="keep all sequences [default]")
        parser.add_argument('-V', '--version', action='version', version=program_version_message)

        # Process arguments
        args = parser.parse_args()

        bam_fp = args.bam
        primer_fp = args.primers
        out_fp = args.out
        wiggle = args.wiggle
        maskbam = args.maskbam
        primeronly = args.primeronly        

        logfile = "primer_masking.log"
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s %(levelname)-8s %(message)s',
                            datefmt='%m/%d/%Y %H:%M:%S',
                            filename=logfile,
                            filemode='w')

        output = "primer_masking.tsv"

        samdata = pysam.AlignmentFile(bam_fp.name, "rb")

        if not out_fp:
            out_fp = "%s_primerMasked.bam" % (os.path.splitext(os.path.basename(samdata.filename.decode("utf-8")))[0])
            
        samout = _primer_mask(samdata, primer_fp, wiggle, maskbam, primeronly, out_fp)

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
        profile_filename = 'asap.maskPrimers_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    main()

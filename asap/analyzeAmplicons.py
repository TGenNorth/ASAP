#!/usr/bin/env python3
# encoding: utf-8
'''
asap.analyzeAmplicons -- Align and interpret amplicon sequencing reads

asap.analyzeAmplicons

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
import skbio
import pkg_resources

from asap import dispatcher
from asap import assayInfo
from asap import __version__
from asap import cmdParser


DEBUG = 1
TESTRUN = 0
PROFILE = 0

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
        if not isinstance(argv, argparse.Namespace):
            sys.argv.extend(argv)
            pass

#     program_name = os.path.basename(sys.argv[0])
#     program_version = "v%s" % __version__
#     program_build_date = str(__updated__)
#     program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
#     if __name__ == '__main__':
#         program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
#     else:
#         program_shortdesc = __doc__.split("\n")[1]
#     #program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
#     program_license = '''%s
#
#   Created by TGen North on %s.
#   Copyright 2015 TGen North. All rights reserved.
#
#   Available for academic and research use only under a license
#   from The Translational Genomics Research Institute (TGen)
#   that is free for non-commercial use.
#
#   Distributed on an "AS IS" basis without warranties
#   or conditions of any kind, either express or implied.
#
#   ***You are running the most current development version,
#   use at your own risk***
#
# USAGE
# ''' % (program_shortdesc, str(__date__))

    try:
        # Setup argument parser
        # parser = argparse.ArgumentParser(prog='Tgen-ASAP', description=program_license, formatter_class=argparse.RawDescriptionHelpFormatter)
        # subparsers = parser.add_subparsers(help='sub-command help')
        # parser_analyzeAmplicons = subparsers.add_parser('analyzeAmplicons', help='analyzeAmplicons help')
        # required_group = parser_analyzeAmplicons.add_argument_group("required arguments")
        # required_group.add_argument("-n", "--name", required=True, help="name for this run. [REQUIRED]")
        # required_group.add_argument("-j", "--json", metavar="FILE", required=True, type=argparse.FileType('r'), help="JSON file of assay descriptions. [REQUIRED]")
        # reads_bams_group = required_group.add_mutually_exclusive_group(required=True)
        # reads_bams_group.add_argument("-r", "--read-dir", dest="rdir", metavar="DIR", help="directory of read files to analyze.")
        # reads_bams_group.add_argument("--bam-dir", dest="bdir", metavar="DIR", help="directory of bam files to analyze.")
        # optional_group = parser_analyzeAmplicons.add_argument_group("optional arguments")
        # optional_group.add_argument("-mf", "--memory-free-buffer", dest="mem_free_buffer", default=500, help="amount of memory to keep free when job manager is 'none' or an unknown selection.")
        # optional_group.add_argument("-o", "--out-dir", dest="odir", metavar="DIR", help="directory to write output files to. [default: `pwd`]")
        # optional_group.add_argument("-s", "--submitter", dest="job_manager", default="SLURM", help="cluster job submitter to use (PBS, SLURM, SLURM_NO_ARRAY, SGE, TASK, none). [default: SLURM]")
        # optional_group.add_argument("--submitter-args", dest="sargs", metavar="ARGS", help="additional arguments to pass to the job submitter, enclosed in \"\".")
        # optional_group.add_argument("--submitter-outpath", dest="submitter_logs_outpath", default="./SlurmJobFiles/", help="informs submitter of a preferred place to put logging files.")
        # optional_group.add_argument("--ts-mailto", dest="tsmailto", default="", help="Email to set TS_MAILTO to, this is useful only when using TASK.")
        # optional_group.add_argument("--smor", action="store_true", default=False, help="perform SMOR analysis with overlapping reads. [default: False]")
        # optional_group.add_argument("-w", "--whole-genome", action="store_true", dest="wholegenome", default=False, help="JSON file uses a whole genome reference, so don't write out the consensus, depth, and proportion arrays for each sample")
        # optional_group.add_argument("--allele-output-threshold", dest="allele_min_reads", default=8, type=int, help="cutoff of # of reads below which allels for amino acids and nucleotide alleles will not be output [default: 8]")
        # optional_group.add_argument("--remove-dups", action="store_true", dest="remove_dups", default=False, help="remove duplicate reads (probably only makes sense for WGS or WMTS data. [default: False]")
        # optional_group.add_argument("--subsample", dest="numreads", type=int, help="Subsample all read files/pairs down to NUMREADS reads/pairs.")
        # optional_group.add_argument("--min-base-qual", dest="bqual", default=5, type=int, help="what is the minimum base quality score (BQS) to use a position (Phred scale, i.e. 10=90, 20=99, 30=99.9 accuracy")
        # if min-base-qual=0 and primer-mask=False, then this should behave almost identically to before my edits, there may be some 'negligable' differences
        # trim_group = parser_analyzeAmplicons.add_argument_group("read trimming options")
        # on_off_group = trim_group.add_mutually_exclusive_group()
        # on_off_group.add_argument("--trim", default="bbduk", help="perform adapter trimming on reads. [default: bbduk]. NOTE: cannot trim-primers if not using bbduk")
        # on_off_group.add_argument("--no-trim", dest="trim", action="store_false", help="do not perform adapter trimming.")
        # trim_group.add_argument("--adapter-sequences", dest="adapters", default=pkg_resources.resource_filename(__name__,'illumina_adapters_all.fasta'), help="location of the adapter sequence file to use for trimming. [default: <ASAP install dir>/asap/illumina_adapters_all.fasta]")
        # trim_group.add_argument("--trim-primers", dest="primers", default=False, help="location of primer file to use for primer trimming. NOTE: Not possible if not using bbduk")
        # trim_group.add_argument("-q", "--qual", nargs="?", const="SLIDINGWINDOW:5:20", help="perform quality trimming [default: False], optional parameter can be used to customize quality trimming parameters to trimmomatic. [default: SLIDINGWINDOW:5:20]")
        # trim_group.add_argument("-l", "--minlen", metavar="LEN", default=80, type=int, help="minimum read length to keep after trimming. [default: 80]")
        # align_group = parser_analyzeAmplicons.add_argument_group("read mapping options")
        # align_group.add_argument("-a", "--aligner", default="bowtie2", help="aligner to use for read mapping, supports bowtie2, novoalign, and bwa. [default: bowtie2]")
        # align_group.add_argument("--aligner-args", dest="aargs", metavar="ARGS", default='', help="additional arguments to pass to the aligner, enclosed in \"\".")
        # align_group.add_argument("-d", "--depth", default=100, type=int, help="minimum read depth required to consider a position covered. [default: 100]")
        # align_group.add_argument("-b", "--breadth", default=0.8, type=float, help="minimum breadth of coverage required to consider an amplicon as present. [default: 0.8]")
        # align_group.add_argument("-p", "--proportion", default=-1, type=float, help="minimum proportion required to call a mutation at a given locus. [default: 0.1]") #Set default to -1 so I can tell if user actually set the value or not
        # align_group.add_argument("-m", "--mutation-depth", dest="mutdepth", default=5, type=int, help="minimum number of reads required to call a mutation at a given locus. [default: 5]")
        # align_group.add_argument("-i", "--identity", dest="percid", default=0, type=float, help="minimum percent identity required to align a read to a reference amplicon sequence. [default: 0]")
        # align_group.add_argument("-il", "--identitylist", dest="percidlist", nargs='+', help="amplicon specific percent identies. Space seperated in amplicon then percentage order.")
        # mask_group = parser_analyzeAmplicons.add_argument_group("primer masking options")
        # mask_group.add_argument("--primer-mask", dest="pmask", default=False, help="location of primer file to use for primer masking.")
        # mask_group.add_argument("--primer-wiggle", dest="wiggle", default=9, type=int, help="how many nucleotides outside the primer window should be used to identify primer sequences")
        # mask_group.add_argument("--primer-mask-bam", dest="pmaskbam", default=True, help="should primer sequences in the alignement file be changed to N for easy viewing (otherwise only base qual is set to 0)")
        # # mask_group.add_argument("--primer-only-bam", dest="ponlybam", default=False, help="Should only sequences with primers be considered when calling variants, currently not implemented")
        # consensus_group = parser_analyzeAmplicons.add_argument_group("consensus calling options")
        # consensus_group.add_argument("--consensus-proportion", default=0.8, type=float, help="minimum proportion required to call at base at that position, else 'N'. [default: 0.8]")
        # consensus_group.add_argument("--fill-gaps", nargs="?", const="n", dest="gap_char", help="fill no coverage gaps in the consensus sequence [default: False], optional parameter is the character to use for filling [defaut: n]")
        # consensus_group.add_argument("--mark-deletions", nargs="?", const="_", dest="del_char", help="fill deletions in the consensus sequence [default: False], optional parameter is the character to use for filling [defaut: _]")
        # parser.add_argument("-V", "--version", action="version", version=program_version_message)
        # parser.add_argument("-D", "--debug", action="store_true", default=False, help="turn on debugging mode")

        # Process arguments
        # print(type(argv))
        if isinstance(argv, argparse.Namespace):
            args = argv
            pass
        else:
            args = cmdParser.parser.parse_args(argv)

        run_name = args.name
        json_filename = dispatcher.expandPath(args.json.name)
        read_dir = args.rdir
        bam_dir = args.bdir
        out_dir = args.odir
        trim = args.trim
        if not args.primers:
            primer_seqs = args.primers
        else:
            primer_seqs = dispatcher.expandPath(args.primers)
        primer_mask_file = args.pmask
        if primer_mask_file:
            primer_mask_file = dispatcher.expandPath(args.pmask)
        wiggle = args.wiggle
        pmaskbam = args.pmaskbam
        ponlybam = False # args.ponlybam
        base_qual = args.bqual
        qual = args.qual
        minlen = args.minlen
        aligner = args.aligner
        aligner_args = args.aargs
        depth = args.depth
        breadth = args.breadth
        proportion = args.proportion
        percid = args.percid
        percidlist = args.percidlist
        mutdepth = args.mutdepth
        adapters = dispatcher.expandPath(args.adapters)
        dispatcher.job_manager = args.job_manager.upper()
        dispatcher.job_manager_args = args.sargs
        dispatcher.submitter_logs_outpath = args.submitter_logs_outpath
        dispatcher.mem_free_buffer = args.mem_free_buffer
        dispatcher.mail_to = args.tsmailto
        smor = args.smor
        allele_min_reads = args.allele_min_reads
        debug = args.debug
        wholegenome = args.wholegenome
        remove_dups = args.remove_dups
        subsample = args.numreads
        con_prop = args.consensus_proportion
        fill_gap_char = args.gap_char
        fill_del_char = args.del_char

        if primer_seqs != False and trim != "bbduk":
            response = input(
                "\nPrimer trimming requested using a trimmer other than bbduk, this is not supported functionality.\nContinue using %s as a trimmer without primer trimming [N]? " % trim)
            if not re.match('^[Yy]',response):
                print("Operation cancelled!")
                quit()

        if not out_dir:
            out_dir = os.getcwd()
        if not (read_dir or bam_dir):
            read_dir = os.getcwd()

        out_dir = dispatcher.expandPath(out_dir)
        if read_dir:
            read_dir = dispatcher.expandPath(read_dir)
        if bam_dir:
            bam_dir = dispatcher.expandPath(bam_dir)

        if os.path.exists(out_dir):
            response = input(
                "\nOutput folder %s already exists!\nFiles in it may be overwritten!\nShould we continue anyway [N]? " % out_dir)
            if not re.match('^[Yy]', response):
                print("Operation cancelled!")
                quit()
        else:
            os.makedirs(out_dir)

        logfile = os.path.join(out_dir, "asap.log")
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s %(levelname)-8s %(message)s',
                            datefmt='%m/%d/%Y %H:%M:%S',
                            filename=logfile,
                            filemode='w')

        dispatcher.slurm_outpath = os.path.join(out_dir, dispatcher.slurm_outpath)
        logging.debug("dispatcher.slurm_outpath = %s" % dispatcher.slurm_outpath)
        if os.path.exists(dispatcher.slurm_outpath):
            for oldJobFile in os.listdir(dispatcher.slurm_outpath):
                os.remove(os.path.join(dispatcher.slurm_outpath, oldJobFile))

        logging.info("Combining reads in %s and JSON file: %s for run: %s. Trim=%s Qual=%s" % (read_dir, json_filename, run_name, trim, qual))
        assay_list = assayInfo.parseJSON(args.json.name)

        bam_list = []
        output_files = []
        final_jobs = []
        xml_dir = os.path.join(out_dir, "xml")
        if not os.path.exists(xml_dir):
            os.makedirs(xml_dir)

        if bam_dir:
            bam_list = dispatcher.findBams(bam_dir)

        if read_dir:
            #reference = assayInfo.generateReference(assay_list)
            ref_fasta = os.path.join(out_dir, "reference.fasta")
            skbio.io.registry.write(assayInfo.generateReference(assay_list), 'fasta', ref_fasta)

            index_job = dispatcher.indexFasta(ref_fasta, aligner)
            read_list = dispatcher.findReads(read_dir)
            for read in read_list:
                if (not read.reads):
                    #TODO: write out appropriate xml for samples with empty read files so they show up in results
                    continue
                if trim != False: #if trimming has not been turned off by --no-trim
                    trimmed_reads = dispatcher.trimAdapters(*read, outdir=out_dir, adapters=adapters, quality=qual, minlen=minlen, dependency=index_job, trimmer=trim, primers=primer_seqs)
                    if subsample:
                        subsampled_reads = dispatcher.subsampleReads(trimmed_reads.sample, trimmed_reads.reads, outdir=out_dir, number=subsample, dependency=trimmed_reads.jobid)
                        (bam_file, job_id) = dispatcher.alignReadsToReference(subsampled_reads.sample, subsampled_reads.reads, ref_fasta, out_dir, jobid=subsampled_reads.jobid, aligner=aligner, args=aligner_args, remove_dups=remove_dups)
                    else:
                        (bam_file, job_id) = dispatcher.alignReadsToReference(trimmed_reads.sample, trimmed_reads.reads, ref_fasta, out_dir, jobid=trimmed_reads.jobid, aligner=aligner, args=aligner_args, remove_dups=remove_dups)
                else: #if trimming has been turned off go straight to aligning
                    if subsample:
                        subsampled_reads = dispatcher.subsampleReads(*read, outdir=out_dir, number=subsample, dependency=index_job)
                        (bam_file, job_id) = dispatcher.alignReadsToReference(subsampled_reads.sample, subsampled_reads.reads, ref_fasta, out_dir, jobid=subsampled_reads.jobid, aligner=aligner, args=aligner_args, remove_dups=remove_dups)
                    else:
                        (bam_file, job_id) = dispatcher.alignReadsToReference(read.sample, read.reads, ref_fasta, out_dir, jobid=index_job, aligner=aligner, args=aligner_args, remove_dups=remove_dups)
                bam_list.append((read.sample, bam_file, job_id))
            if trim == "bbduk":
                if not os.path.isdir(out_dir+"/trimmed/STATS/"):
                    os.makedirs(out_dir+"/trimmed/STATS/")
                    pass
                f = open(out_dir+"/trimmed/STATS/info.txt", "w+")
                f.write("If run with paired reads the stats information for each pair is in [read_pair_name]*R1*STATS, if run with single reads the stats information for the single read is in [read_name]*R1*STATS. _primers if is the stats for trimming primers, no _primers if stats for trimming adapters")
                f.close()
        for sample, bam, job in bam_list:
            (xml_file, job_id) = dispatcher.processBam(sample, json_filename, bam, xml_dir, job, depth, breadth, proportion, percid, mutdepth, primer_mask_file, wiggle, pmaskbam, ponlybam, base_qual, smor, wholegenome, debug, allele_min_reads, con_prop, fill_gap_char, fill_del_char, percidlist)

            output_files.append(xml_file)
            final_jobs.append(job_id)

        (final_output, job) = dispatcher.combineOutputFiles(run_name, xml_dir, out_dir, final_jobs)
        print("All jobs are submitted, the final job id is: %s. Output will be in %s when ready." % (job, final_output))

        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception as e:
        if DEBUG or TESTRUN:
            raise(e)
        indent = len(cmdParser.program_name) * " "
        sys.stderr.write(cmdParser.program_name + ": " + repr(e) + "\n")
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

#!/usr/bin/env python3
# encoding: utf-8
'''
asap.cmdParser -- argparse setup for asap

asap.cmdParser

@author:     Darrin Lemmer & Caleb Johnson

@copyright:  2015,2019 TGen North. All rights reserved.

@license:    ACADEMIC AND RESEARCH LICENSE -- see ../LICENSE

@contact:    dlemmer@tgen.org
'''

import argparse
from asap import assayInfo
from asap import __version__
import os
import sys
import re
import logging
import skbio
import pkg_resources

from asap import analyzeAmplicons
from asap import bamProcessor
from asap import outputCombiner
from asap import formatOutput
from asap import prepareJSONInput
from asap import reformatXML

__all__ = []
__updated__ = '2020-07-20'
__date__ = '2015-06-04'

program_name = "asap"
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

parser = argparse.ArgumentParser(prog=program_name, description=program_license, formatter_class=argparse.RawDescriptionHelpFormatter)

def main(argv=None):

        subparsers = parser.add_subparsers(help='Sub-commands list, used via '+program_name+' <sub-command>')

        # Analyze Amplicons
        parser_analyzeAmplicons = subparsers.add_parser('analyzeAmplicons', help=program_name+' analyzeAmplicons --help')
        required_group = parser_analyzeAmplicons.add_argument_group("required arguments")
        required_group.add_argument("-n", "--name", required=True, help="name for this run. [REQUIRED]")
        required_group.add_argument("-j", "--json", metavar="FILE", required=True, type=argparse.FileType('r'), help="JSON file of assay descriptions. [REQUIRED]")
        reads_bams_group = required_group.add_mutually_exclusive_group(required=True)
        reads_bams_group.add_argument("-r", "--read-dir", dest="rdir", metavar="DIR", help="directory of read files to analyze.")
        reads_bams_group.add_argument("--bam-dir", dest="bdir", metavar="DIR", help="directory of bam files to analyze.")
        optional_group = parser_analyzeAmplicons.add_argument_group("optional arguments")
        optional_group.add_argument("-mf", "--memory-free-buffer", dest="mem_free_buffer", default=500, help="amount of memory to keep free when job manager is 'none' or an unknown selection.")
        optional_group.add_argument("-o", "--out-dir", dest="odir", metavar="DIR", help="directory to write output files to. [default: `pwd`]")
        optional_group.add_argument("-s", "--submitter", dest="job_manager", default="SLURM", help="cluster job submitter to use (PBS, SLURM, SLURM_NO_ARRAY, SGE, TASK, none). [default: SLURM]")
        optional_group.add_argument("--submitter-args", dest="sargs", metavar="ARGS", help="additional arguments to pass to the job submitter, enclosed in \"\".")
        optional_group.add_argument("--submitter-outpath", dest="submitter_logs_outpath", default="./SlurmJobFiles/", help="informs submitter of a preferred place to put logging files.")
        optional_group.add_argument("--ts-mailto", dest="tsmailto", default="", help="Email to set TS_MAILTO to, this is useful only when using TASK.")
        optional_group.add_argument("--smor", action="store_true", default=False, help="perform SMOR analysis with overlapping reads. [default: False]")
        optional_group.add_argument("-w", "--whole-genome", action="store_true", dest="wholegenome", default=False, help="JSON file uses a whole genome reference, so don't write out the consensus, depth, and proportion arrays for each sample")
        optional_group.add_argument("--allele-output-threshold", dest="allele_min_reads", default=8, type=int, help="cutoff of # of reads below which allels for amino acids and nucleotide alleles will not be output [default: 8]")
        optional_group.add_argument("--remove-dups", action="store_true", dest="remove_dups", default=False, help="remove duplicate reads (ony makes sense for WGS or WMTS data. [default: False]")
        trim_group = parser_analyzeAmplicons.add_argument_group("read trimming options")
        on_off_group = trim_group.add_mutually_exclusive_group()
        on_off_group.add_argument("--trim", default="bbduk", help="perform adapter trimming on reads. [default: bbduk]. NOTE: cannot trim-primers if not using bbduk")
        on_off_group.add_argument("--no-trim", dest="trim", action="store_false", help="do not perform adapter trimming.")
        trim_group.add_argument("--adapter-sequences", dest="adapters", default=pkg_resources.resource_filename(__name__,'illumina_adapters_all.fasta'), help="location of the adapter sequence file to use for trimming. [default: <ASAP install dir>/asap/illumina_adapters_all.fasta]")
        trim_group.add_argument("--trim-primers", dest="primers", default=False, help="location of primer file to use for primer trimming. NOTE: Not possible if not using bbduk")
        trim_group.add_argument("-q", "--qual", nargs="?", const="SLIDINGWINDOW:5:20", help="perform quality trimming [default: False], optional parameter can be used to customize quality trimming parameters to trimmomatic. [default: SLIDINGWINDOW:5:20]")
        trim_group.add_argument("-l", "--minlen", metavar="LEN", default=80, type=int, help="minimum read length to keep after trimming. [default: 80]")
        align_group = parser_analyzeAmplicons.add_argument_group("read mapping options")
        align_group.add_argument("-a", "--aligner", default="bowtie2", help="aligner to use for read mapping, supports bowtie2, novoalign, and bwa. [default: bowtie2]")
        align_group.add_argument("--aligner-args", dest="aargs", metavar="ARGS", default='', help="additional arguments to pass to the aligner, enclosed in \"\".")
        align_group.add_argument("-d", "--depth", default=100, type=int, help="minimum read depth required to consider a position covered. [default: 100]")
        align_group.add_argument("-b", "--breadth", default=0.8, type=float, help="minimum breadth of coverage required to consider an amplicon as present. [default: 0.8]")
        align_group.add_argument("-p", "--proportion", type=float, help="minimum proportion required to call a mutation at a given locus. [default: 0.1]") #Don't explicitly set default because I need to be certain whether user set the value
        align_group.add_argument("-m", "--mutation-depth", dest="mutdepth", default=5, type=int, help="minimum number of reads required to call a mutation at a given locus. [default: 5]")
        align_group.add_argument("-i", "--identity", dest="percid", default=0, type=float, help="minimum percent identity required to align a read to a reference amplicon sequence. [default: 0]")
        align_group.add_argument("-il", "--identitylist", dest="percidlist", nargs='+', help="amplicon specific percent identies. Space seperated in amplicon then percentage order.")
        mask_group = parser_analyzeAmplicons.add_argument_group("primer masking options")
        mask_group.add_argument("--primer-mask", dest="pmask", default=False, help="location of primer file to use for primer masking.")
        mask_group.add_argument("--primer-wiggle", dest="wiggle", default=9, type=int, help="how many nucleotides outside the primer window should be used to identify primer sequences")
        mask_group.add_argument("--primer-mask-bam", dest="pmaskbam", default=True, help="should primer sequences in the alignement file be changed to N for easy viewing (otherwise only base qual is set to 0)")
        # mask_group.add_argument("--primer-only-bam", dest="ponlybam", default=False, help="Should only sequences with primers be considered when calling variants, currently not implemented")
        optional_group.add_argument("--min_base_qual", dest="bqual", default=5, type=int, help="what is the minimum base quality score (BQS) to use a position (Phred scale, i.e. 10=90, 20=99, 30=99.9 accuracy")
        # if min_base_qual=0 and primer-mask=False, then this should behave almost identically to before my edits, there may be some 'negligable' differences
        consensus_group = parser_analyzeAmplicons.add_argument_group("consensus calling options")
        consensus_group.add_argument("--consensus-proportion", default=0.8, type=float, help="minimum proportion required to call at base at that position, else 'N'. [default: 0.8]")
        consensus_group.add_argument("--fill-gaps", nargs="?", const="n", dest="gap_char", help="fill no coverage gaps in the consensus sequence [default: False], optional parameter is the character to use for filling [defaut: n]")
        consensus_group.add_argument("--mark-deletions", nargs="?", const="_", dest="del_char", help="fill deletions in the consensus sequence [default: False], optional parameter is the character to use for filling [defaut: _]")
        parser_analyzeAmplicons.add_argument("-V", "--version", action="version", version=program_version_message)
        parser_analyzeAmplicons.add_argument("-D", "--debug", action="store_true", default=False, help="turn on debugging mode")
        parser_analyzeAmplicons.set_defaults(func=analyzeAmplicons.main)

        # Bam Processor
        parser_bamProcessor = subparsers.add_parser('bamProcessor', help=program_name+' bamProcessor --help')
        required_group = parser_bamProcessor.add_argument_group("required arguments")
        required_group.add_argument("-j", "--json", metavar="FILE", required=True, type=argparse.FileType('r'), help="JSON file of assay descriptions. [REQUIRED]")
        required_group.add_argument("-b", "--bam", metavar="FILE", required=True, type=argparse.FileType('rb'), default=sys.stdin, help="BAM file to analyze. [REQUIRED]")
        #required_group.add_argument("-r", "--ref", metavar="FILE", required=True, help="reference fasta file, should already be indexed. [REQUIRED]")
        #parser.add_argument("-o", "--out-dir", dest="odir", metavar="DIR", help="directory to write output files to. [default: `pwd`]")
        # TODO: (argparse file type and optional. default to stdout)
        #parser.add_argument("-n", "--name", help="sample name, if not provided it will be derived from BAM file")
        parser_bamProcessor.add_argument("-d", "--depth", default=100, type=int, help="minimum read depth required to consider a position covered. [default: 100]")
        parser_bamProcessor.add_argument("--breadth", default=0.8, type=float, help="minimum breadth of coverage required to consider an amplicon as present. [default: 0.8]")
        parser_bamProcessor.add_argument("-p", "--proportion", default=0.1, type=float, help="minimum proportion required to call a mutation at a given locus. [default: 0.1]")
        parser_bamProcessor.add_argument("-m", "--mutation-depth", dest="mutdepth", default=5, type=int, help="minimum number of reads required to call a mutation at a given locus. [default: 5]")
        identity_group = parser_bamProcessor.add_argument_group("identity filter options")
        identity_group.add_argument("-i", "--identity", dest="percid", default=0, type=float, help="minimum percent identity required to align a read to a reference amplicon sequence. [default: 0]")
        identity_group.add_argument("-il", "--identitylist", dest="percidlist", nargs='+', help="amplicon specific percent identies. Space seperated in amplicon then percentage order.")
        keep_discarded_group = identity_group.add_mutually_exclusive_group()
        keep_discarded_group.add_argument("-k", "--keep", action="store_true", default=False, help="keep filtered reads. [default: True]")
        keep_discarded_group.add_argument("--no-keep", action="store_false", dest="keep", help="discard filtered reads.")
        merge_group = identity_group.add_mutually_exclusive_group()
        merge_group.add_argument("--merge", action="store_true", default=False, help="merge paired reads. [default: False]")
        merge_group.add_argument("--no-merge", action="store_false", dest="merge", help="do not merge paired reads.")
        parser_bamProcessor.add_argument("-s", "--smor", action="store_true", default=False, help="perform SMOR analysis with overlapping reads. [default: False]")
        parser_bamProcessor.add_argument("-V", "--version", action="version", version=program_version_message)
        parser_bamProcessor.add_argument("-D", "--debug", action="store_true", default=False, help="write <sample_name>.log file with debugging information")
        parser_bamProcessor.add_argument("-w", "--whole-genome", action="store_true", dest="wholegenome", default=False, help="JSON file uses a whole genome reference, so don't write out the consensus, depth, and proportion arrays for each sample")
        parser_bamProcessor.add_argument("--allele-output-threshold", dest="allele_min_reads", default=8, type=int, help="cutoff of # of reads below which allels for amino acids and nucleotide alleles will not be output [default: 8]")
        parser_bamProcessor.add_argument('-o', '--out', metavar="FILE", type=argparse.FileType('w'), default=sys.stdout, help="output filename [default: stdout]")
        parser_bamProcessor.add_argument("--output-format", type=str.lower, choices=('xml', 'json'), default='xml', help="output format [default: xml]")
        parser_bamProcessor.add_argument("--primer-mask", dest = "pmask", default=False, help="Location of primer file to use for primer masking.")
        parser_bamProcessor.add_argument("--primer-wiggle", dest="wiggle", default=9, type=int, help="How many nucleotides outside the primer window should be used to identify primer sequences")
        parser_bamProcessor.add_argument("--primer-mask-bam", dest="pmaskbam", default=True, help="Should primer sequences in the alignement file be changed to N")
        parser_bamProcessor.add_argument("--primer-only-bam", dest="ponlybam", default=True, help="Should only sequences with primers be considered when calling variants")
        parser_bamProcessor.add_argument("--min_base_qual", dest="bqual", default=5, type=int, help="What is the minimum base quality score to use a position (phred scale, i.e. 10=90, 20=99, 30=99.9, accuraccy")
        parser_bamProcessor.add_argument("--consensus-proportion", default=0.8, type=float, help="minimum proportion required to call at base at that position, else 'N'. [default: 0.8]")
        parser_bamProcessor.add_argument("--fill-gaps", nargs="?", const="n", dest="gap_char", help="fill no coverage gaps in the consensus sequence [default: False], optional parameter is the character to use for filling [defaut: n]")
        parser_bamProcessor.add_argument("--mark-deletions", nargs="?", const="_", dest="del_char", help="fill deletions in the consensus sequence [default: False], optional parameter is the character to use for filling [defaut: _]")
        parser_bamProcessor.set_defaults(func=bamProcessor.main)

        #Output Combiner
        parser_outputCombiner = subparsers.add_parser('outputCombiner', help=program_name+' outputCombiner --help')
        required_group = parser_outputCombiner.add_argument_group("required arguments")
        required_group.add_argument("-n", "--name", required=True, help="name for this run. [REQUIRED]")
        required_group.add_argument("-x", "--xml-dir", dest="xdir", metavar="DIR", required=True, help="directory containing XML files to combine. [REQUIRED]")
        parser_outputCombiner.add_argument("-o", "--out", metavar="FILE", help="file to write final output to. [default: ./{name}_analysis.xml]")
        parser_outputCombiner.add_argument('-V', '--version', action='version', version=program_version_message)
        parser_outputCombiner.set_defaults(func=outputCombiner.main)

        # Format Output
        parser_formatOutput = subparsers.add_parser('formatOutput', help=program_name+' formatOutput --help')
        required_group = parser_formatOutput.add_argument_group("required arguments")
        required_group.add_argument("-s", "--stylesheet", metavar="FILE", required=True, help="XSLT stylesheet to use for transforming the output. [REQUIRED]")
        required_group.add_argument("-x", "--xml", metavar="FILE", required=True, help="XML output file to transform. [REQUIRED]")
        parser_formatOutput.add_argument("-o", "--out", dest="out", metavar="FILE", help="output file to write.")
        parser_formatOutput.add_argument("-d", "--outdir", dest="out_dir", metavar="DIR", help="output directory to write files to.")
        parser_formatOutput.add_argument("-t", "--text", action="store_true", default=False, help="output plain text.")
        parser_formatOutput.add_argument('-V', '--version', action='version', version=program_version_message)
        parser_formatOutput.set_defaults(func=formatOutput.main)

        #prepareJSONInput
        parser_prepareJSONInput = subparsers.add_parser('prepareJSONInput', help=program_name+' prepareJSONInput --help')
        required_group = parser_prepareJSONInput.add_argument_group("required arguments")
        exclusive_group = required_group.add_mutually_exclusive_group(required=True)
        exclusive_group.add_argument("-f", "--fasta", metavar="FILE", help="fasta file containing amplicon sequences.")
        exclusive_group.add_argument("-x", "--excel", metavar="FILE", help="Excel file of assay data.")
        required_group.add_argument("-o", "--out", metavar="FILE", required=True, help="output JSON file to write. [REQUIRED]")
        parser_prepareJSONInput.add_argument("-w", "--worksheet", help="Excel worksheet to use, the first one in the file will be used if not specified")
        parser_prepareJSONInput.add_argument('-v', '--version', action='version', version=program_version_message)
        parser_prepareJSONInput.set_defaults(func=prepareJSONInput.main)

        #reformatXML
        parser_reformatXML = subparsers.add_parser('reformatXML', help=program_name+' reformatXML --help')
        parser_reformatXML.add_argument("-x", "--xml", required=True, help="ASAP output XML file to reformat. [REQUIRED]")
        parser_reformatXML.add_argument("-V", "--version", action="version", version=program_version_message)
        parser_reformatXML.set_defaults(func=reformatXML.main)


        args = parser.parse_args()

        try:
            args.func(args)
            pass
        except AttributeError as e:
            print("\nFor help use: "+program_name+" -h")


if __name__ == ' __main__':
    sys.exit(main())

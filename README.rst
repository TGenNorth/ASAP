.. |copy|   unicode:: U+000A9 .. COPYRIGHT SIGN

Amplicon Sequencing Analysis Pipeline (ASAP)
========================================

OVERVIEW:
---------
The Amplicon Sequencing Analysis Pipeline (ASAP) is a highly customizable, automated way to examine amplicon sequencing data. The important details of the amplicon targets are described in a text-based input file written in JavaScript Object Notation (JSON) [1]_. This data includes the target name, amplicon sequence (or sequences in the case of gene variant assays), any known SNPs or regions of interest (ROIs) within the target, and what the presence of this target or SNP signifies. This file can be hand-generated or created from an Excel spreadsheet using a provided template and Python script. The sequenced reads are processed by performing adapter, and optionally, quality trimming using Trimmomatic [2]_, and then aligned to the reference amplicon sequences extracted from the JSON file using one of several alignment packages (BWA-MEM [3]_, bowtie2 [4]_, and NovoAlign [5]_ are currently supported). The resulting BAM [6]_ files are analyzed with a custom Python script using the pysam [7]_ and scikit-bio [8]_ libraries to aid in analysis. This script combines the alignment data in the BAM file with the assay data in the JSON file and interprets the results. The output is an XML file with complete details for each assay against each sample. These details include number of reads aligning to each target, any SNPs found above a user-defined threshold, and the nucleotide distribution at each of these SNP positions. For ROI assays, the output includes the sequence distribution at each of the regions of interest -- both the DNA sequences and translated into amino acid sequences. Also, each assay target is assigned a significance if it meets the requirements laid out in the JSON file (i.e. a particular SNP or amino acid change is present) To make this output easier for the user to interpret, a number of XSLT [9]_ stylesheets are provided for transforming the XML output into other, more readable formats, including Excel spreadsheets, web pages, and PDF documents. Additionally, the use of XSLT stylesheets allows for multiple different views of the same data, from clinical summaries showing only the most important or relevant results to full researcher summaries containing all of the data.

USAGE:
------
**1) Generating JSON File**

Can be generated from Excel spreadsheet template, or for simple cases, directly from multifasta file.

typical usage: ``prepareJSONInput -x <EXCEL_FILE> -o <OUTPUT_JSON_FILE>``

full usage: ``prepareJSONInput [-h] (-f FILE | -x FILE) -o FILE [-w WORKSHEET] [-V]``

asap.prepareJSONInput -- Create a JSON input file for ASAP from a multifasta or Excel spreadsheet

optional arguments:
  -h, --help            show this help message and exit
  -w WORKSHEET, --worksheet WORKSHEET
                        Excel worksheet to use, the first one in the file will
                        be used if not specified
  -V, --version         show program's version number and exit

required arguments:
  -f FILE, --fasta FILE
                        fasta file containing amplicon sequences.
  -x FILE, --excel FILE
                        Excel file of assay data.
  -o FILE, --out FILE   output JSON file to write. [REQUIRED]


**2) Running ASAP**

typical usage: ``analyzeAmplicons -n <RUN_NAME> -j <PATH_TO_JSON_FILE> -r <DIRECTORY_OF_READ_FILES> -o <OUTPUT_DIRECTORY> <other options>``

``<RUN_NAME>`` can be whatever you want, the final output file will be: ``<OUTPUT_DIRECTORY>/<RUN_NAME>_analysis.xml``

You can also change the depth (default 100), proportion (default 0.1), breadth (default 0.8) filters using the ``-d``, ``-p`` and ``-b`` options

full usage: ``analyzeAmplicons [-h] -n NAME -j JSON [-r DIR | --bam-dir DIR] [-o DIR] [-s JOB_MANAGER] [--submitter-args ARGS] [--smor] [--trim | --no-trim] [-s ADAPTERS] [-q [QUAL]] [-m LEN] [-a ALIGNER] [--aligner-args ARGS] [-d DEPTH] [--breadth BREADTH] [-p PROPORTION] [-i PERCID] [-V]``

asap.analyzeAmplicons -- Align and interpret amplicon sequencing reads

optional arguments:
  -h, --help            show this help message and exit
  -V, --version         show program's version number and exit

required arguments:
  -n NAME, --name NAME  name for this run. [REQUIRED]
  -j JSON, --json JSON  JSON file of assay descriptions. [REQUIRED]

optional arguments:
  -r DIR, --read-dir DIR
                        directory of read files to analyze.
  --bam-dir DIR         directory of bam files to analyze.
  -o DIR, --out-dir DIR
                        directory to write output files to. [default: `pwd`]
  -s JOB_MANAGER, --submitter JOB_MANAGER
                        cluster job submitter to use (PBS, SLURM, SGE, none).
                        [default: PBS]
  --submitter-args ARGS
                        additional arguments to pass to the job submitter,
                        enclosed in "".
  --smor                perform SMOR analysis with overlapping reads.
                        [default: False]

read trimming options:
  --trim                perform adapter trimming on reads. [default: True]
  --no-trim             do not perform adapter trimming.
  -s ADAPTERS, --adapter-sequences ADAPTERS
                        location of the adapter sequence file to use for
                        trimming.
  -q QUAL, --qual QUAL
                        perform quality trimming [default: False], optional
                        parameter can be used to customize quality trimming
                        parameters to trimmomatic. [default:
                        SLIDINGWINDOW:5:20]
  -m LEN, --minlen LEN  minimum read length to keep after trimming. [default:
                        80]

read mapping options:
  -a ALIGNER, --aligner ALIGNER
                        aligner to use for read mapping, supports bowtie2,
                        novoalign, and bwa. [default: bowtie2]
  --aligner-args ARGS   additional arguments to pass to the aligner, enclosed
                        in "".
  -d DEPTH, --depth DEPTH
                        minimum read depth required to consider a position
                        covered. [default: 100]
  -b BREADTH, --breadth BREADTH     
                        minimum breadth of coverage required to consider an
                        amplicon as present. [default: 0.8]
  -p PROPORTION, --proportion PROPORTION
                        minimum proportion required to call a SNP at a given
                        position. [default: 0.1]
  -i PERCID, --identity PERCID
                        minimum percent identity required to align a read to a
                        reference amplicon sequence. [default: 0]

This command will ultimately generate the xml file. To convert this into more better things, run:


**3) Formatting Output**

typical usage ``formatOutput -s <XSLT_FILE> -x <XML_OUTPUT_FILE> -o <MAIN_OUTPUT_FILE_TO_WRITE>``

This will generate all the html files, which you can open directly in your web browser. Some xslt files are available in the ``output_transforms`` directory.

full usage: ``formatOutput [-h] -s FILE -x FILE [-o FILE] [-t] [-V]``

asap.formatOutput -- Apply an XSLT transformation on the XML output to generate a more user-friendly output

optional arguments:
  -h, --help            show this help message and exit
  -t, --text            output plain text
  -V, --version         show program's version number and exit

required arguments:
  -s FILE, --stylesheet FILE
                        XSLT stylesheet to use for transforming the output.
                        [REQUIRED]
  -x FILE, --xml FILE   XML output file to transform. [REQUIRED]
  -o FILE, --out FILE   output file to write. [REQUIRED]


INSTALLATION:
-------------

See the included "INSTALL" document.

DEPENDENCIES:
-------------

For information about external tools that are required, or can be
utilized, and those versions that have been tested to work with ASAP,
refer to the included "INSTALL" document.

LICENSE:
--------

Copyright |copy| The Translational Genomics Research Institute See the
included "LICENSE" document.

CONTACT:
--------

Darrin Lemmer (dlemmer@tgen.org)
| TGen North
| 3051 W Shamrell Blvd Ste 106
| Flagstaff, AZ 86001-9435

REFERENCES:
-----------

.. [1] JSON: http://www.ecma-international.org/publications/files/ECMA-ST/ECMA-404.pdf
.. [2] Trimmomatic: Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.
.. [3] BWA-MEM: http://bio-bwa.sourceforge.net - There’s a publication for BWA-SW, and BWA short read aligner, but not for BWA-MEM. Maybe the short read aligner paper should be referenced here? The details are at this link.
.. [4] Bowtie2: Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012, 9:357-359.
.. [5] NovoAlign: http://www.novocraft.com - seems there should be a better reference, but I haven’t found one.
.. [6] SAM format/SAMtools: Li, Heng et al. “The Sequence Alignment/Map Format and SAMtools.” Bioinformatics 25.16 (2009): 2078–2079. PMC. Web. 9 Nov. 2015.
.. [7] Pysam: https://github.com/pysam-developers/pysam
.. [8] Scikit-bio: http://scikit-bio.org
.. [9] XSLT: http://www.w3.org/TR/xslt20/

#!/usr/bin/env python3
# encoding: utf-8
'''
asap.reformatXML -- re-write the XML output file to be sorted by assay instead of sample

asap.reformatXML

@author:     Jonathon Todd

@copyright:  2018,2019 TGen North. All rights reserved.

@license:    ACADEMIC AND RESEARCH LICENSE -- see ../LICENSE

@contact:    dlemmer@tgen.org
'''

import xml.etree.ElementTree as ET
import copy
import sys
import os
import argparse

__all__ = []
__version__ = 0.6
__date__ = '2018-04-11'
__updated__ = '2019-01-15'

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

def removeChildren (node):
    victims = []
    for child in node:
        victims.append(child)
    for child in victims:
        node.remove(child)


def main(argv=None): # IGNORE:C0111
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        if not isinstance(argv, argparse.Namespace):
            sys.argv.extend(argv)
            pass

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
        # parser = argparse.ArgumentParser(description=program_license, formatter_class=argparse.RawDescriptionHelpFormatter)
        # parser.add_argument("-x", "--xml", required=True, help="ASAP output XML file to reformat. [REQUIRED]")
        # parser.add_argument("-V", "--version", action="version", version=program_version_message)

        # Process arguments
        if isinstance(argv, argparse.Namespace):
            args = argv
            pass
        else:
            args = asapParser.parser.parse_args(argv)

        tree = ET.parse(args.xml)
        root = tree.getroot()
        #print root[0][0].attrib["name"]

        tree2 = copy.deepcopy(tree)
        #tree2.getroot().clear()
        #tree2._setroot(tree.getroot())
        root2 = tree2.getroot()
        #victims = []
        #for child in root2:
         #   victims.append(child)
        #for child in victims:
         #   root2.remove(child)
        removeChildren(root2)

        # insert all assays from first sample
        sample1 = root.find('sample')
        #print sample1
        for assay in sample1:
            tmp = copy.deepcopy(assay)
            removeChildren(tmp)
            root2.append(tmp)

        #insert all sample into those assays
        # get all samples
        samples = root.findall('sample')
        assays2 = root2.findall('assay')
        #print assays2
        for assay in assays2:
            for sample in samples:
                xpath = "assay[@name='%s']" % assay.attrib["name"]
                if sample.find(xpath) != None:
                    tmp = copy.deepcopy(sample)
                    removeChildren(tmp)
                    xpath2 = "sample[@bam_file='%s']/assay[@name='%s']" % (sample.attrib["bam_file"], assay.attrib["name"])
                    origAssays = root.findall(xpath2)
                    #print origAssays
                    for a in origAssays:
                        for child in a:
                            #print "a child"
                            tmp.append(child)
                            assay.append(tmp)


        #this section edits the names of all unknown snps to be Unknown:{POSITION}{CALL}
        #TODO: Let's make this optional in the future, this is a hack for a specific problem, and we do not want to always do this
        for assay in assays2:
            if assay.attrib["type"] == 'SNP' or assay.attrib["type"] == 'mixed':
                snps = assay.findall('.//snp')
                for  snp in snps:
                    if snp.attrib["name"] == "unknown":
                        snp.set("name","unknown:" + snp.attrib["position"] + snp.find('snp_call').text)

        tree2.write(args.xml.split(".")[0]+"_Reformated.xml")

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
        profile_filename = 'asap.reformatXML_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())

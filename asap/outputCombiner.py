#!/usr/bin/env python3
# encoding: utf-8
'''
asap.outputCombiner -- Combine multiple sample-level XML output files into one run-level XML output

asap.outputCombiner 

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
from xml.etree import ElementTree

from asap import dispatcher
from asap import __version__ 

__all__ = []
__date__ = '2015-07-29'
__updated__ = '2019-01-15'

DEBUG = 1
TESTRUN = 0
PROFILE = 0

def _write_xml(root, xml_file):
    from xml.dom import minidom
    dom = minidom.parseString(ElementTree.tostring(root))
    output = open(xml_file, 'w')
    output.write('\n'.join([line for line in dom.toprettyxml(indent=' '*2).split('\n') if line.strip()]))
    output.close()
    return xml_file

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
        parser = argparse.ArgumentParser(description=program_license, formatter_class=argparse.RawDescriptionHelpFormatter)
        required_group = parser.add_argument_group("required arguments")
        required_group.add_argument("-n", "--name", required=True, help="name for this run. [REQUIRED]")
        required_group.add_argument("-x", "--xml-dir", dest="xdir", metavar="DIR", required=True, help="directory containing XML files to combine. [REQUIRED]")
        parser.add_argument("-o", "--out", metavar="FILE", help="file to write final output to. [default: ./{name}_analysis.xml]")
        parser.add_argument('-V', '--version', action='version', version=program_version_message)

        # Process arguments
        args = parser.parse_args()

        run_name = args.name
        xml_dir = args.xdir
        out_file = args.out

        xml_dir = dispatcher.expandPath(xml_dir)

        root_node = ElementTree.Element("analysis", {'run_name':run_name})
        
        sorted_files = sorted(os.listdir(xml_dir))    
        for file in sorted_files:
            is_xml = re.search('(.*)\.xml$', file, re.IGNORECASE)
            if is_xml:
                sample_tree = ElementTree.parse(os.path.join(xml_dir, file))
                sample_node = sample_tree.getroot()
                root_node.append(sample_node)

        if not out_file:
            out_file = os.path.join(".", run_name+"_analysis.xml")
            
        _write_xml(root_node, out_file)

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
        profile_filename = 'asap.outputCombiner_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())

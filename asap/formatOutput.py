#!/usr/bin/env python3
# encoding: utf-8
'''
asap.formatOutput -- Apply an XSLT transformation on the XML output to generate a more user-friendly output

asap.formatOutput

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
from xml.etree import ElementTree
import lxml.etree as ET

import asap.dispatcher as dispatcher

__all__ = []
__version__ = 0.1
__date__ = '2015-07-29'
__updated__ = '2015-07-29'

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
        required_group.add_argument("-s", "--stylesheet", metavar="FILE", required=True, help="XSLT stylesheet to use for transforming the output. [REQUIRED]")
        required_group.add_argument("-x", "--xml", metavar="FILE", required=True, help="XML output file to transform. [REQUIRED]")
        required_group.add_argument("-o", "--out", dest="out", metavar="FILE", help="output file to write. [REQUIRED]")
        parser.add_argument('-V', '--version', action='version', version=program_version_message)

        # Process arguments
        args = parser.parse_args()

        stylesheet = args.stylesheet
        xml_file = args.xml
        out_file = args.out
        run_name = ""

        match = re.search('^(.*)_analysis.xml$', xml_file)
        if match:
            run_name = match.group(1)
            
        if run_name and not os.path.exists(run_name):
            os.makedirs(run_name)

        dom = ET.parse(xml_file)
        xslt = ET.parse(stylesheet)
        transform = ET.XSLT(xslt)
        newdom = transform(dom)
        
        output = open(out_file, 'wb')
        output.write(ET.tostring(newdom, pretty_print=True, xml_declaration=True, encoding='UTF-8'))
        output.close()

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

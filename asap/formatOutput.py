#!/usr/bin/env python3
# encoding: utf-8
'''
asap.formatOutput -- Apply an XSLT transformation on the XML output to generate a more user-friendly output

asap.formatOutput

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
import lxml.etree as ET

from asap import dispatcher
from asap import __version__

__all__ = []
__date__ = '2015-07-29'
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

def distinct_values(context, values):
    return list(set(values))

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
        parser.add_argument("-o", "--out", dest="out", metavar="FILE", help="output file to write.")
        parser.add_argument("-d", "--outdir", dest="out_dir", metavar="DIR", help="output directory to write files to.")
        parser.add_argument("-t", "--text", action="store_true", default=False, help="output plain text.")
        parser.add_argument('-V', '--version', action='version', version=program_version_message)

        # Process arguments
        args = parser.parse_args()

        stylesheet = args.stylesheet
        xml_file = args.xml
        out_file = args.out
        out_dir = args.out_dir
        text = args.text

        if not out_dir:
            match = re.search('^(.*)_analysis.xml$', xml_file)
            if match:
                out_dir = match.group(1)
            
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir)

        ns = ET.FunctionNamespace("http://pathogen.tgen.org/ASAP/functions")
        ns['distinct-values'] = distinct_values

        #parser = ET.XMLParser(huge_tree=True)
        #dom = ET.parse(xml_file, parser=parser)
        dom = ET.parse(xml_file)
        xslt = ET.parse(stylesheet)
        transform = ET.XSLT(xslt)
        newdom = transform(dom)
        
        if out_file:
            output = open(out_file, 'wb')
            if text:
                output.write(newdom)
            else:
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

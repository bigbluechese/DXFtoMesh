#!/usr/bin/python
'''
Code Written by Jeff Peterson and Kerry Wang (2016)
Takes a DXF file and converts it into segments (which can be arcs or lines) and
vertexes that can then be converted into a CrysMAS or Cats2D mesh. 

Requirements:
- Python 2.7x
- dxfgrabber 0.7.5
- matplotlib (if plotting)
- DXFtoSegments.py
- HelperFunctions.py

Requirements for Testing:
- DXFTestSuite.py (if using this for testing)
- DXFTests directory of DXF test files

Usage:
This script takes input from the command line and the usage for which can be
found by running this script with the '-h' flag afterwards specificaly,
'./DXFtoSegments.py -h'. 
'''

import DXFTestSuite
from DXFtoSegments import DXFGeometry
import CreateStandards
import unittest
import argparse
from matplotlib import pyplot as plt
import sys

def test_suite(verbose=False, dxf_test=True):
    '''Runs the test suite'''
    if verbose:
        verbosity = 2
    else:
        verbosity = 1
    suites = []
    suites.append(unittest.TestLoader().loadTestsFromTestCase(DXFTestSuite.TestVertex))
    suites.append(unittest.TestLoader().loadTestsFromTestCase(DXFTestSuite.TestVertexList))
    if dxf_test:
        suites.append(unittest.TestLoader().loadTestsFromTestCase(DXFTestSuite.TestDXFGeometry))
        suites.append(unittest.TestLoader().loadTestsFromTestCase(DXFTestSuite.DXFTestCases))
        suites.append(unittest.TestLoader().loadTestsFromTestCase(DXFTestSuite.DXFtoCats2DTests))
        dxf1 = DXFGeometry('./DXFTests/DXFTest1.dxf', testing=True)
        dxf1.display()
        dxf2 = DXFGeometry('./DXFTests/DXFTest2.dxf', testing=True)
        dxf2.display()
        dxf3 = DXFGeometry('./DXFTests/DXFTest3.dxf', testing=True)
        dxf3.display()
        dxf4 = DXFGeometry('./DXFTests/DXFTest4.dxf', testing=True)
        dxf4.display()
    print
    alltests = unittest.TestSuite(suites)
    unittest.TextTestRunner(verbosity=verbosity).run(alltests)
    # Show the plotted DXF file if available
    plt.show()

# NOTE: Argument parsing requires Python 2.7x or higher
# Parses command-line input arguemtns
help_string = '''Reads a DXF file and then converts it to a suitable form
                 for use in creating computational meshes for CrysMAS and 
                 Cats2D'''
parser = argparse.ArgumentParser(description=help_string)
# Specify the DXF file (required)
help_string = '''Specify the DXF file to convert or specify \'test\' to run
                 the test suite'''
parser.add_argument('dxf_file', action='store', type=str, 
                    metavar='dxf_file or \'test\'', help=help_string)
# Create verbose mode
parser.add_argument('-v', '--verbose', action='store_true')
# Specify whether information should be output to CrysMAS and optionally specify a new filename
help_string = '''Turn the DXF file into a CrysMAS .pcs mesh file. A file 
                 name for the .pcs file can optionally be specified'''
parser.add_argument('-c', '--crysmas', nargs='?', metavar='file',
                    const=True, help=help_string)
# Skip DXF tests if option is passed
help_string = '''Skips the DXF tests if testing mode is activated'''
parser.add_argument('--nodxf', action='store_true', help=help_string)

# Create new standard DXF files for tests
help_string = '''Creates new standards for tests to compare to when in test mode. Use with caution'''
parser.add_argument('--newstandard', action='store_true', help=help_string)

args = parser.parse_args()

# Specify testing mode from the command line
if args.dxf_file == 'test':
    # Create new standards if specified
    if args.newstandard:
        print 'Creating new standards for tests'
        CreateStandards.save_standards(num_tests=4)
    else:
        test_suite(verbose=args.verbose, dxf_test=not(args.nodxf))
    # Exit immediately after testing is complete
    sys.exit()
else: #Otherwise create a DXF geometry object
    dxf = DXFGeometry(args.dxf_file, verbose=args.verbose)

if args.crysmas:
    print 'Make a crysmas geometry'
    if args.crysmas != True:
        fname = args.crysmas
    
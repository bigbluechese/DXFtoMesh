#!/usr/bin/env python
'''
Code Written by Jeff Peterson and Kerry Wang (2016)
Takes a DXF file and converts it into segments (which can be arcs or lines) and
vertexes that can then be converted into a CrysMAS or Cats2D mesh. 

Requirements:
- Python 2.7x
- dxfgrabber 0.7.5 (0.8.0 is incompatible without patch)
- SciPy Stack
- sortedcontainers 1.5.3
- DXFtoSegments.py
- HelperFunctions.py

Requirements for Testing:
- DXFTestSuite.py
- CreateStandards.py
- DXFTests directory of DXF test files

Usage:
This script takes input from the command line and the usage for which can be
found by running this script with the '-h' flag afterwards specificaly,
'./DXFtoSegments.py -h'. 
'''

from __future__ import print_function
import DXFTestSuite
from DXFtoSegments import DXFGeometry
import CreateStandards
import unittest
import argparse
import matplotlib as mpl
mpl.use('TkAgg')
from matplotlib import pyplot as plt
import os
import sys
import tarfile
import c2d_premesh_v5
import pexpect_c2dmesh_v2


def test_suite(verbose=False, dxf_test=True):
    '''Runs the test suite'''
    if verbose:
        verbosity = 2
    else:
        verbosity = 1
    suites = []
    suites.append(unittest.TestLoader().loadTestsFromTestCase(DXFTestSuite.TestVertex))
    suites.append(unittest.TestLoader().loadTestsFromTestCase(DXFTestSuite.TestVertexList))
    suites.append(unittest.TestLoader().loadTestsFromTestCase(DXFTestSuite.LineIntersectTests))
    if dxf_test:
        suites.append(unittest.TestLoader().loadTestsFromTestCase(DXFTestSuite.TestDXFGeometry))
        suites.append(unittest.TestLoader().loadTestsFromTestCase(DXFTestSuite.DXFTestCases))
        suites.append(unittest.TestLoader().loadTestsFromTestCase(DXFTestSuite.DXFtoCats2DTests))
        suites.append(unittest.TestLoader().loadTestsFromTestCase(DXFTestSuite.DXFtoCrysMASTests))
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

def make_dist():
    '''Creates a tar file for the distribution of this code'''
    # Find out where this code is being run from
    root_path = os.path.dirname(os.getcwd())
    base_path = os.path.basename(os.getcwd())
    os.chdir('..')
    # Set file types to ignore
    ignore_types = ['.pyc', '.pcs', '.dwl', '.dwl2']
    tar_name = base_path+'.tar.gz'
    with tarfile.open(tar_name, 'w:gz') as tar:
        for root, dirs, files in os.walk(base_path):
            # Ignore hidden files and .pyc files
            files = [f for f in files if (not f[0] == '.') and (not [ext for ext in ignore_types if ext in f[-4:]])]
            # Ignore hidden directories
            dirs[:] = [d for d in dirs if not d[0] == '.']
            # Add these files to the tar file
            for f in files:
                print('adding', os.path.join(root, f))
                tar.add(os.path.join(root, f))
        print('packaging code as {} at {}'.format(tar_name, root_path))
        tar.close()
    os.chdir('./DXFtoMesh')

# NOTE: Argument parsing requires Python 2.7x or higher
# Parses command-line input arguemtns
help_string = '''Reads a DXF file and then converts it to a suitable form for
                 use in creating computational meshes for CrysMAS and Cats2D'''
parser = argparse.ArgumentParser(description=help_string)

# Specify the DXF file (required)
help_string = '''Specify the DXF file to convert, specify \'test\' to run the
                 test suite, or \'dist\' to create the distribution of this
                 program'''
parser.add_argument('dxf_file', action='store', type=str, 
                    metavar='dxf_file or \'test\'', help=help_string)

# Create verbose mode
parser.add_argument('-v', '--verbose', action='store_true')

# Specify whether information should be output to CrysMAS and optionally specify
# a new filename
help_string = '''Turn the DXF file into a CrysMAS .pcs mesh file. A file name
                 for the .pcs file can optionally be specified'''
parser.add_argument('--crysmas', nargs='?', metavar='file',
                    const=True, help=help_string)

# Specify whether information should be output to Cats2D and what the length
#  scale should be for non-dimensionalization
help_string = '''Turn the DXF file into a Cats2D mesh. The length scale for
                 non-dimensionalization must also be specified'''
parser.add_argument('--cats2d', action='store', metavar='len_scale', 
                    help=help_string)

# Specify the execution name for Cats2D
help_string = '''Optionally set the path for running Cats2D'''
parser.add_argument('--c2dpath', action='store', default='cats2d.x',
                    help=help_string)

# Plot Cats2D regions
help_string = '''Optionally plot Cats2D regions as they are found. The delay
                 between plotting events can optionally be specified'''
parser.add_argument('--c2dplot', action='store', nargs='?', metavar='delay',
                    const=0.2, help=help_string)

# Specify the units of the DXF file
help_string = 'Specify the units for the DXF file'
parser.add_argument('--units', action='store', default='mm', help=help_string)

# Pickle the DXFGeometry object
help_string = 'Create a pickle file of the DXFGeoemtry object'
parser.add_argument('--pickle', nargs='?', const=True, help=help_string)

# Skip DXF tests if option is passed
help_string = '''Skips the DXF tests if testing mode is activated'''
parser.add_argument('--nodxf', action='store_true', help=help_string)

# Create new standard DXF files for tests
help_string = '''Creates new standards for tests to compare to when in test
                mode. Use with caution'''
parser.add_argument('--newstandard', action='store_true', help=help_string)

args = parser.parse_args()

# Specify testing mode from the command line
if args.dxf_file == 'test':
    # Create new standards if specified
    if args.newstandard:
        print('Creating new standards for tests')
        CreateStandards.save_standards(args, num_tests=4)
    else:
        test_suite(verbose=args.verbose, dxf_test=not(args.nodxf))
    # Exit immediately after testing is complete
    sys.exit()
# Create distribution from the command line
elif args.dxf_file == 'dist':
    make_dist()
# Otherwise create a DXF geometry object
else: 
    dxf = DXFGeometry(args.dxf_file, verbose=args.verbose)

# Create a CrysMAS file
if args.crysmas:
    if args.crysmas != True:
        crys_file = args.crysmas
    else:
        crys_file = None
    dxf.output_to_crysmas(dxf_units=args.units, f_name=crys_file)

# Create a Cats2D mesh
if args.cats2d:
    if args.c2dplot:
        plotting_args = {'plotting':True, 'plotting_pause':args.c2dplot}
    else:
        plotting_args = {}
    vertex_list,edge_list,bulge_list = dxf.cats2d_convert(len_scale=6)
    mesh = c2d_premesh_v5.C2DMesh(vertex_list, edge_list, **plotting_args)
    pexpect_c2dmesh_v2.make_c2d_mesh(mesh, args.c2dpath, working_dir=dxf.work_dir)

# Create a pickle file
if args.pickle:
    if args.pickle != True:
        dxf.pickle(f_name=args.pickle)
    else:
        dxf.pickle()

    
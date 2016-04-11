#!/usr/bin/python
'''
Function for creating new test standards for the DXF-specific tests. Create new
standards with caution as overwriting the correct standard may hide errors that
the testing suite won't identify.
'''

import numpy as np
import DXFtoSegments
from matplotlib import pyplot as plt
import pickle

def save_standards(num_tests=4):
    '''Creates standards from DXF test files to be used for testing purposes.'''
    geometries = {}
    # Open the DXF files and create DXFGeometry objects
    for i in range(num_tests):
        directory = './DXFTests/'
        file_name = 'DXFTest{}'.format(i+1)
        extension = '.dxf'
        path = directory+file_name+extension
        dxf = DXFtoSegments.DXFGeometry(path, verbose=False, testing=True)
        geometries[file_name] = dxf

        # Pickle the segments set
        ppath = '{}DXFTest{}_segments.set'.format(directory, i+1)
        f = open(ppath, 'wb')
        pickle.dump(dxf.segments, f)
        f.close()

        # Create Cats2D-compatible information
        if i == 3:
            verts, edges, bulges = dxf.cats2d_convert()
            cats_ppath = '{}DXFTest{}_cats2d.pick'.format(directory, i+1)
            f_cats = open(cats_ppath, 'wb')
            pickle.dump((verts, edges, bulges), f_cats)
            f_cats.close()
            print '\t\tcats2d_convert data stored for {}'.format(file_name)

        # Now check whether it was successful
        f_new = open(ppath, 'rb')
        file_segs = pickle.load(f_new)
        if file_segs == dxf.segments:
            print '\tData stored for {} in file {}'.format(file_name, ppath)
        else:
            raise Exception('Something went wrong when pickling the data')
        f_new.close()

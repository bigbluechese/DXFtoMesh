#!/usr/bin/python

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
        dxf = DXFtoSegments.DXFGeometry(path, verbose=False)
        geometries[file_name] = dxf

        # Pickle the segments set
        ppath = './DXFTests/DXFTest{}_segments.set'.format(i+1)
        f = open(ppath, 'wb')
        pickle.dump(dxf.segments, f)
        f.close()

        # Now check whether it was successful
        f_new = open(ppath, 'rb')
        file_segs = pickle.load(f_new)
        if file_segs == dxf.segments:
            print 'Data stored for {} in file {}'.format(file_name, ppath)
        else:
            raise Exception('Something went wrong when pickling the data')
        f_new.close()
#!/usr/bin/env python
'''
Example for the usage of the DXFGeometry object within python
'''
import numpy as np
import DXFtoSegments
from matplotlib import pyplot as plt

tests = 4

dxf_objects = []

for i in range(tests):
    file_name = './DXFTests/DXFTest{}.dxf'.format(i+1)
    dxf_objects.append(DXFtoSegments.DXFGeometry(file_name))

dxf1 = dxf_objects[0]
dxf2 = dxf_objects[1]
dxf3 = dxf_objects[2]
dxf4 = dxf_objects[3]

clamshell = DXFtoSegments.DXFGeometry('./DXFTests/DXFTest_Clamshellv5.dxf')

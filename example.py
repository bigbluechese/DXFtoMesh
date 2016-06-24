#!/usr/bin/env python
'''
Example for the usage of the DXFGeometry object within python
'''
from __future__ import print_function
import numpy as np
import DXFtoSegments
import lineintersect
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

ITest2 = DXFtoSegments.DXFGeometry('./DXFTests/IntersectTest2.dxf')
intersections = lineintersect.find_intersections(ITest2.verts.vertices, verbose=True, tol=0.3)

# ITest2.display()
# plt.show()

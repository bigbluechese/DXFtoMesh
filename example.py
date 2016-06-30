#!/usr/bin/env python
'''
Example for the usage of the DXFGeometry object within python
'''
from __future__ import print_function
import numpy as np
import DXFtoSegments
import lineintersect
from matplotlib import pyplot as plt
import time


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
repeats = 10
t0 = time.time()
for i in range(repeats):
    intersections1 = lineintersect.find_intersections_brute(clamshell.verts.vertices, verbose=False)
t1 = time.time()
for i in range(repeats):
    intersections2 = lineintersect.find_intersections(clamshell.verts.vertices, verbose=False)
t2 = time.time()

msg = '\n\nTIME FOR BRUTE FORCE: {}\nTIME FOR SWEEP LINE: {}'
print(msg.format((t1-t0)/repeats, (t2-t1)/repeats))
print('Number of intersections:')
print('Brute force: {} \t\t Sweep line: {}'.format(len(intersections1), len(intersections2)))

# ITest2.display()
# plt.show()

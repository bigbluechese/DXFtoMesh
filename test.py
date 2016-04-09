import numpy as np
import DXFtoSegments
from matplotlib import pyplot as plt

tests = 4

dxf_objects = []

for i in range(tests):
    file_name = './DXFTests/DXFTest{}.dxf'.format(i+1)
    dxf_objects.append(DXFtoSegments.DXFGeometry(file_name, verbose=True))

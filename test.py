import numpy as np
import DXFtoSegments
from matplotlib import pyplot as plt

dxf1 = DXFtoSegments.DXFGeometry('./DXFTests/DXFTest1.dxf', verbose=True)
dxf2 = DXFtoSegments.DXFGeometry('./DXFTests/DXFTest2.dxf')
dxf3 = DXFtoSegments.DXFGeometry('./DXFTests/DXFTest3.dxf')
dxf4 = DXFtoSegments.DXFGeometry('./DXFTests/DXFTest3.dxf')

point1 = list(dxf1.verts.coordinates)[6]
point2 = eval(list(dxf1.verts.verticies[str(point1)].connected)[0])

for i, s in enumerate(dxf1.segments):
    if point1 in s[0] and point2 in s[0]:
        ind = i
        coords = s[0]

seg = list(dxf1.segments)[ind]
info = seg[1]
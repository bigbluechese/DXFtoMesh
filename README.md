# DXFtoMesh
Written by Jeff Peterson and Kerry Wang (2016)

This program is used to read a DXF file and convert the entities in the DXF file
into forms that can be used to construct either CrysMAS or Cats2D meshes.
Specifically, this program takes a DXF file, breaks it up into its entities, and
then breaks those entities into "segments" that are described by the coordinates
of the endpoints and then various pieces of information abou the curvature of
the segment.

Example usage:

python MeshGenerator.py ./DXFTests/DXFTest_Clamshellv5.dxf --crysmas

REQUIREMENTS:
- Python 2.7x
- dxfgrabber 0.7.5
- matplotlib (if plotting)

REQUIRED FILES:
- DXFtoSegments.py
- HelperFunctions.py

FILES REQUIRED FOR TESTING:
- DXFTestSuite.py
- CreateStandards.py
- ./DXFTests/ (directory)

For use inside of other scripts, these modules may be used by importing the
DXFtoSegments module. For example,

import DXFtoSegments
dxf_geometry = DXFtoSegments.DXFGeometry(file_name)

will create a DXFGeometry object which can then output to CrysMAS by invoking
the output_to_crysmas() method.
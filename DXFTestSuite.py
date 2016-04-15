#!/usr/bin/python

'''
Code Written by Jeff Peterson (2016)
A suite of tests for checking the functionality of the DXFtoSegments module.
'''

import unittest
from DXFtoSegments import VertexList, Vertex, DXFGeometry
import dxfgrabber
import math
import numpy as np
from HelperFunctions import bulge_to_arc, approx, ccw_angle_diff
import pickle
import filecmp

class TestVertex(unittest.TestCase):
    '''
    Test case for the VertexList class. Note that the self.v instance should be
    separate for each test case so they shouldn't influence each other.
    '''
    def setUp(self):
        x = -15.9384
        y = 1.2345e+02
        self.v = Vertex((x,y))

    def tearDown(self):
        self.v = None

    def test_integer_creation(self):
        '''Create vertex using integers'''
        x = 1
        y = -2
        dummy_v = Vertex((x,y))
        self.assertEqual(dummy_v.x, x, 'incorrect x-value stored')
        self.assertEqual(dummy_v.y, y, 'incorrect y-value stored')
        self.assertEqual(dummy_v.id, str((float(x),float(y))),
                            'vertex id is not string of tuple')

    def test_decimal_creation(self):
        '''Create vertex using decimals'''
        x = 0.1
        y = -0.2
        dummy_v = Vertex((x,y))
        self.assertEqual(dummy_v.x, x, 'incorrect x-value stored')
        self.assertEqual(dummy_v.y, y, 'incorrect y-value stored')
        self.assertEqual(dummy_v.id, str((float(x),float(y))),
                            'vertex id is not string of tuple')

    def test_scientifc_creation(self):
        '''Create vertex using scientific numbers'''
        x = 1.0e+01
        y = -2.4e1
        dummy_v = Vertex((x,y))
        self.assertEqual(dummy_v.x, x, 'incorrect x-value stored')
        self.assertEqual(dummy_v.y, y, 'incorrect y-value stored')
        self.assertEqual(dummy_v.id, str((float(x),float(y))),
                            'vertex id is not string of tuple')

    def test_connect(self):
        '''Add connection information to vertex'''
        new_vertex = (0.,1.)
        self.v.con(str(new_vertex))
        check = str(new_vertex) in self.v.connected
        self.assertTrue(check, 'new vertex not connected')

    def test_connect_multiple(self):
        '''Add multiple connections to vertex'''
        new_vertex1 = (0.,1.)
        new_vertex2 = (-2., -2.)
        self.v.con(str(new_vertex1))
        self.v.con(str(new_vertex2))
        check = (str(new_vertex1) in self.v.connected) and \
                    (str(new_vertex2) in self.v.connected)
        self.assertTrue(check, 'new verticies not connected')

    def test_disconnect(self):
        '''Remove connection information from vertex'''
        new_vertex = (0,1) #First connect a vertex
        self.v.con(str(new_vertex))
        self.v.discon(str(new_vertex)) #Now disconnect it
        check = new_vertex in self.v.connected
        self.assertFalse(check, 'new vertex still connected')


class TestVertexList(unittest.TestCase):
    '''
    \tTest case for the VertexList class
    '''
    def setUp(self):
        self.v_grid = set([]) #Create an empty set of vertexes to be tested
        x_steps, y_steps = 10, 10 #Range/domain of verticies
        self.x_low, self.x_high, self.y_low, self.y_high = 0., 1., 0., 1.
        self.v_list = VertexList()
        for x in np.linspace(self.x_low, self.x_high, x_steps):
            for y in np.linspace(self.y_low, self.y_high, y_steps):
                self.v_grid.add((x, y)) #Add each vertex

    def tearDown(self):
        self.v_grid = None
        self.v_list = None

    def test_add_verticies(self):
        '''Add new vetecies'''
        for v in self.v_grid:
            self.v_list.add(v)
        check = self.v_list.coordinates == self.v_grid
        self.assertTrue(check, 'vertex list set is not equal to set to be added')

    def test_connect_verticies(self):
        '''Connect verticies'''
        # First add all verticies
        for v in self.v_grid:
            self.v_list.add(v)
        # Now connect corners together
        corners = set([])
        for x in np.linspace(self.x_low, self.x_high, 2):
            for y in np.linspace(self.y_low, self.y_high, 2):
                corners.add(str((x,y)))
        # Connect vertex to every other vertex except itself
        for c1 in corners:
            for c2 in corners:
                if c1 != c2:
                    self.v_list.connect(c1, c2)
                else:
                    # Check that a vertex cannot be connected to itself
                    with self.assertRaises(RuntimeError):
                        self.v_list.connect(c1, c1)
        # Check connectivity
        checks = []
        for c in corners:
            other_corners = corners.copy()
            other_corners.remove(c)
            check = self.v_list.verticies[str(c)].connected == other_corners
            if not check:
                checks.append(check)

        msg = '{} of {} verticies did not have proper connectivity'.format(len(checks),
                                                                        len(corners))
        self.assertFalse(bool(checks), msg)


    def test_disconnect_verticies(self):
        '''Disconnect verticies'''
        # First add all verticies
        for v in self.v_grid:
            self.v_list.add(v)
        # Now all connect corners together
        corners = set([])
        for x in np.linspace(self.x_low, self.x_high, 2):
            for y in np.linspace(self.y_low, self.y_high, 2):
                corners.add(str((x,y)))
        # Connect vertex to every other vertex except itself
        for c1 in corners:
            for c2 in corners:
                if c1 != c2:
                    self.v_list.connect(c1, c2)
                else:
                    pass
        # Now disconnect the an arbitrary corner from another arbitrary corner
        new_connections = corners.copy()
        corner1 = new_connections.pop()
        corner2 = new_connections.pop()
        self.v_list.disconnect(corner1, corner2)
        check = (self.v_list.verticies[corner1].connected == new_connections) \
                and (self.v_list.verticies[corner2].connected == new_connections)
        self.assertTrue(check)

    def test_remove_vertex(self):
        '''Remove vertex from list without connections'''
        # First add them
        for v in self.v_grid:
            self.v_list.add(v)
        # Create a copy and remove a vertex
        v_grid_new = self.v_grid.copy()
        v_remove = v_grid_new.pop()
        self.v_list.remove(v_remove)
        check = (self.v_list.coordinates == v_grid_new) \
                 and (v_remove not in self.v_list.coordinates)
        self.assertTrue(check, 'vertex was not removed properly')
        # Now check if the vertex has been removed from the dict
        with self.assertRaises(KeyError):
            self.v_list.verticies[str(v_remove)]

    def test_move_vertex(self):
        '''Move a vertex'''
        # First add all verticies
        for v in self.v_grid:
            self.v_list.add(v)
        # Now move the lowest vertex to a negative value
        old_position = (self.x_low, self.y_low)
        new_position = (self.x_low - (self.x_high - self.x_low), 
                        self.y_low - (self.y_high - self.y_low))
        self.v_list.move_vertex(old_position, new_position)
        check = (new_position in self.v_list.coordinates) \
                and (old_position not in self.v_list.coordinates)
        self.assertTrue(check, 'vertex was not moved')

class TestDXFGeometry(unittest.TestCase):
    '''\tTest case for the DXFGeometry class'''
    def setUp(self):
        # Open DXF file for all tests
        self.test_dxf = dxfgrabber.readfile('./DXFTests/DXFTest2.dxf')
        self.empty_dxfgeom = DXFGeometry('empty_test', testing=True)

    def tearDown(self):
        # Get rid of DXF file
        self.test_dxf = None
        self.empty_dxfgeom = None

    def check_verticies(self, tup, dxfgeom, entity):
        '''verticies are in vertex list and they are connected'''
        tol = dxfgeom.tol
        start = (approx(tup[0][0], tol), approx(tup[0][1],tol))
        end = (approx(tup[1][0], tol), approx(tup[1][1],tol))
        # Check to make sure start and end are verticies now
        check = (start in dxfgeom.verts.coordinates) and (end in dxfgeom.verts.coordinates)
        msg = '{} and {} were not added to verticies'.format(start, end)
        self.assertTrue(check, msg)
        # Make sure verticies are connected
        check = str(end) in dxfgeom.verts.verticies[str(start)].connected
        msg = '{} not properly connected by line {}'.format(end, entity)
        self.assertTrue(check, msg)

    def test_add_line(self):
        '''add a line to the geometry'''
        # Loop through entities and add only lines
        tol = self.empty_dxfgeom.tol
        for e in self.test_dxf.entities:
            if e.dxftype == 'LINE':
                self.empty_dxfgeom.add_line(e)
                start = (approx(e.start[0], tol=tol), 
                         approx(e.start[1], tol=tol))
                end = (e.end[0], e.end[1])
                self.check_verticies((start, end), self.empty_dxfgeom, e)

    def test_add_arc(self):
        '''add an arc to the geometry'''
        # Loop through entities and add only arcs
        for e in self.test_dxf.entities:
            if e.dxftype == 'ARC':
                self.empty_dxfgeom.add_arc(e)
                tol = self.empty_dxfgeom.tol
                start_angle = math.radians(e.startangle)
                end_angle = math.radians(e.endangle)
                center = (approx(e.center[0], tol=tol), 
                          approx(e.center[1], tol=tol))
                radius = e.radius
                # Calculate bulge and start/stop information
                theta = ccw_angle_diff(start_angle, end_angle)
                bulge = math.tan(theta/4)
                start = (approx(radius*math.cos(start_angle) + center[0], tol=tol), 
                         approx(radius*math.sin(start_angle) + center[1], tol=tol))
                end = (approx(radius*math.cos(end_angle) + center[0], tol=tol), 
                       approx(radius*math.sin(end_angle) + center[1], tol=tol))
                self.check_verticies((start, end), self.empty_dxfgeom, e)

    def test_add_polyline(self):
        '''add a polyline to the geometry'''
        # Loop through entities and add only arcs
        tol = self.empty_dxfgeom.tol
        for e in self.test_dxf.entities:
            if e.dxftype == 'POLYLINE':
                self.empty_dxfgeom.add_polyline(e)
                for i, p in enumerate(e.points):
                    try:
                        p_next = e.points[i+1]
                    except IndexError:
                        if e.is_closed:
                            p_next = e.points[0]
                        else:
                            break
                    start = (approx(p[0], tol=tol),
                             approx(p[1], tol=tol))
                    end = (approx(p_next[0], tol=tol), 
                            approx(p_next[1], tol=tol))
                    self.check_verticies((start, end), self.empty_dxfgeom, e)

    def test_add_entities(self):
        '''add entites from a DXF file'''
        # Add entities
        self.empty_dxfgeom.add_entities(self.test_dxf.entities)
        # There should be 18 lines, arcs, and polylines in the file
        num_ents = 18
        num_added = len(self.empty_dxfgeom.segments)
        check = num_added == num_ents
        msg = '''{} entities were added but {} addable entities exist in test file
              '''.format(num_added, num_ents)
        self.assertTrue(check, msg)

    def test_rem_reversed(self):
        '''remove reversed lines and arcs'''
        # Lines to be added
        lines = set([(( (0,0), (1,0) ), ()),
                     (( (1,0), (0,0) ), ()),
                     (( (0,0), (1,1) ), ())])
        # Arcs to be added
        arcs = set([])
        bulges = [1, -1, 1]
        for l, b in zip(lines, bulges):
            start = l[0][0]
            end = l[0][1]
            segment = bulge_to_arc(start, end, b)
            arcs.add(segment)

        # Lines that are not duplicates:
        true_lines = set([(( (0,0), (1,0) ), ()),
                          (( (0,0), (1,1) ), ())])
        # Arcs that are not duplicates:
        true_arcs = set([bulge_to_arc((0,0), (1,0), -1),
                         bulge_to_arc((0,0), (1,1), 1)])

        # Add all lines and verticies
        self.empty_dxfgeom.segments |= lines | arcs
        for l in lines:
            self.empty_dxfgeom.verts.add(l[0][0])
            self.empty_dxfgeom.verts.add(l[0][1])
            self.empty_dxfgeom.verts.connect(l[0][0], l[0][1])

        # Remove reversed line and arc
        self.empty_dxfgeom.rem_reversed()

        # Check whether it worked
        check = self.empty_dxfgeom.segments == true_lines | true_arcs
        msg = '''Extra line that should have been removed: {}
              '''.format(self.empty_dxfgeom.segments - true_lines | true_arcs)
        self.assertTrue(check, msg)

class DXFTestCases(unittest.TestCase):
    '''Specific test cases for comparing DXF Geometries against known standards'''
    def setUp(self):
        self.test_path = './DXFTests/'
        self.num_tests = 4

    def test_DXFs(self):
        '''DXF geometries are same as saved standards'''
        for i in range(self.num_tests):
            dxf = DXFGeometry('{0}DXFTest{1}.dxf'.format(self.test_path, i+1), testing=True)
            standard = open('{0}DXFTest{1}_segments.set'.format(self.test_path, 
                            i+1), 'rb')
            std_segs = pickle.load(standard)
            check = dxf.segments == std_segs
            msg = '''DXFTest{} Geometry does not match standard geometry'''.format(i)
            self.assertTrue(check, msg)

class DXFtoCats2DTests(unittest.TestCase):
    '''Test cases for converting a DXF geometry into a Cats2D mesh'''
    def setUp(self):
        self.dxf4 = DXFGeometry('./DXFTests/DXFTest4.dxf', testing=True)

    def tearDown(self):
        self.dxf4 = None

    def test_DXFGeomtoCats(self):
        '''Convert DXF geometry for creating Cats2D mesh'''
        standard = open('./DXFTests/DXFTest4_cats2d.pick', 'rb')
        pick_info = pickle.load(standard)
        check = self.dxf4.cats2d_convert(invert_coords=False) == pick_info
        msg = 'Cats2D information for DXFTest4 did not match the saved standard'
        self.assertTrue(check, msg)

class DXFtoCrysMASTests(unittest.TestCase):
    '''Test cases for converting a DXF geometry into a CrysMAS .pcs file'''
    def setUp(self):
        self.dxf = DXFGeometry('./DXFTests/DXFTest_Clamshellv5.dxf')

    def tearDown(self):
        self.dxf = None

    def test_DXFGeomtoCrys(self):
        '''Conver DXF geometry to CrysMAS .pcs file'''
        standard = './DXFTests/DXFTest_Clamshellv5.pcs.const'
        self.dxf.output_to_crysmas(f_name='Testing')
        test_file = './DXFTests/Testing.pcs'
        check = filecmp.cmp(test_file, standard)
        msg = 'File information differs between {} and {}'.format(standard, test_file)
        self.assertTrue(check, msg)

def main():
    print '''Please run the test suite from the MeshGenerator.py file by using
    the command line argument `test` with the optional verbose mode.'''

# Check whether the script is being excuted by itself
if __name__=='__main__':
    main()

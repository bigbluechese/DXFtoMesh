#!/usr/bin/python

'''
Code Written by Jeff Peterson (2016)
A suite of tests for checking the functionality of the DXFtoSegments module.

Required Python Modules:
- DXFtoSegments
- MatPlotLib
'''

import unittest
from DXFtoSegments import VertexList, Vertex, DXFGeometry, drange
import dxfgrabber
import math
#import matplotlib.pyplot as plt

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
        for x in drange(self.x_low, self.x_high, x_steps):
            for y in drange(self.y_low, self.y_high, y_steps):
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
        for x in drange(self.x_low, self.x_high, 2):
            for y in drange(self.y_low, self.y_high, 2):
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
        for x in drange(self.x_low, self.x_high, 2):
            for y in drange(self.y_low, self.y_high, 2):
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
        self.empty_dxfgeom = DXFGeometry('./DXFTests/DXFTest2.dxf', testing=True)

    def tearDown(self):
        # Get rid of DXF file
        self.test_dxf = None
        self.empty_dxfgeom = None

    def test_openDXF(self):
        '''
        Test that opens the file DXFTest2 and creates a DXFGeometry without
        adding any entities
        '''
        # Use testing mode to keep from adding entities
        dxf2 = DXFGeometry('./DXFTests/DXFTest2.dxf', testing=True)
        # Test that the first entity is consistent
        check = dxf2.dxf.entities[0].dxftype == self.test_dxf.entities[0].dxftype
        self.assertTrue(check)
        # Make sure that there are the same number of entities
        check = len(dxf2.dxf.entities) == len(self.test_dxf.entities)
        self.assertTrue(check)

    def check_verticies(self, tup, dxfgeom, entity):
        '''Checks whether verticies are in vertex list and whether they are
        connected'''
        start = tup[0]
        end = tup[1]
        # Check to make sure start and end are verticies now
        check = (start in dxfgeom.verts.coordinates) and (end in dxfgeom.verts.coordinates)
        msg = '{} and {} were not added to verticies'.format(start, end)
        self.assertTrue(check, msg)
        # Make sure verticies are connected
        check = str(end) in dxfgeom.verts.verticies[str(start)].connected
        msg = '{} not properly connected by line {}'.format(end, entity)
        self.assertTrue(check, msg)

    def test_add_line(self):
        '''Tests ability to add a line to the geometry'''
        # Loop through entities and add only lines
        for e in self.test_dxf.entities:
            if e.dxftype == 'LINE':
                self.empty_dxfgeom.add_line(e)
                start = (e.start[0], e.start[1])
                end = (e.end[0], e.end[1])
                self.check_verticies((start, end), self.empty_dxfgeom, e)
                # Make sure lines are added to segments list
                check = ((start, end), ()) in self.empty_dxfgeom.segments
                msg = '{} not in segment list'.format((start, end))
                self.assertTrue(check, msg)
    def test_add_arc(self):
        '''Tests ability to add an arc to the geometry'''
        # Loop through entities and add only arcs
        for e in self.test_dxf.entities:
            if e.dxftype == 'ARC':
                self.empty_dxfgeom.add_arc(e)
                start_angle = math.radians(e.startangle)
                end_angle = math.radians(e.endangle)
                center = (e.center[0], e.center[1])
                radius = e.radius
                # Calculate bulge and start/stop information
                bulge = math.tan((end_angle - start_angle)/4)
                start = (radius*math.cos(start_angle) + center[0], 
                         radius*math.sin(start_angle) + center[1])
                end = (radius*math.cos(end_angle) + center[0], 
                       radius*math.sin(end_angle) + center[1])
                self.check_verticies((start, end), self.empty_dxfgeom, e)
                # Make sure arcs were added to segment list
                check = ((start, end), (bulge, start_angle, end_angle, center, \
                         radius)) in self.empty_dxfgeom.segments
                msg = '{} not in segment list'.format((start, end))
                self.assertTrue(check, msg)
    def test_add_polyline(self):
        '''Tests ability to add an polyline to the geometry'''
        # Loop through entities and add only arcs
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
                    start = (p[0], p[1])
                    end = (p_next[0], p_next[1])
                    self.check_verticies((start, end), self.empty_dxfgeom, e)
                    # Make sure at least straight lines are added to segments
                    # list (check for arcs by inspection)
                    check = ((start, end), ()) in self.empty_dxfgeom.segments
                    msg = '{} not in segment list'.format((start, end))
                    if e.bulge[i] == 0:
                        self.assertTrue(check, msg)

    def test_add_entities(self):
        '''Tests ability to add entites from a DXF file'''
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
        '''Tests the ability to remove reversed lines'''
        # Lines to be added
        lines = set([(( (0,0), (1,0) ), ()),
                     (( (1,0), (0,0) ), ()),
                     (( (0,0), (1,1) ), ())])

        # Lines that are not duplicates
        true_lines = set([(( (0,0), (1,0) ), ()),
                          (( (0,0), (1,1) ), ())])
        # Add all lines and verticies
        self.empty_dxfgeom.segments = lines
        for l in lines:
            self.empty_dxfgeom.verts.add(l[0][0])
            self.empty_dxfgeom.verts.add(l[0][1])
            self.empty_dxfgeom.verts.connect(l[0][0], l[0][1])

        # Remove reversed line
        self.empty_dxfgeom.rem_reversed()

        # Check whether it worked
        check = self.empty_dxfgeom.segments == true_lines
        msg = '''Extra line that should have been removed: {}
              '''.format(self.empty_dxfgeom.segments - true_lines)

    def test_arc_calc(self):
        dxf = DXFGeometry('./DXFTests/DXFTest1.dxf')
        tol = 1e-12 # Tolerance for positional differences
        # Check that segment points are consistent with arc end points
        checks = ''
        for seg in dxf.segments:
            # Skip straight lines
            if not seg[1]:
                continue
            start = seg[0][0]
            end = seg[0][1]
            center = seg[1][3]
            radius = seg[1][4]
            sangle = seg[1][1]
            eangle = seg[1][2]
            calc_start = (radius*math.cos(sangle)+center[0], 
                          radius*math.sin(sangle)+center[1])
            calc_end = (radius*math.cos(eangle)+center[0], 
                        radius*math.sin(eangle)+center[1])
            # Look for differences
            diff = (calc_start[0]-start[0], calc_start[1] - start[1])
            if abs(diff[0]) > tol or abs(diff[1]) > tol:
                msg = '''\tActual starting point: {}
                         Calculated starting point: {}
                         \tDifference: {}
                      '''.format(start, calc_start, diff)
                checks = checks + msg
            diff = (calc_end[0]-end[0], calc_end[1] - end[1])
            tol = 1e-12
            if abs(diff[0]) > tol or abs(diff[1]) > tol:
                msg = '''\tActual end point: {}
                         Calculated end point: {}
                         \tDifference: {}
                      '''.format(end, calc_end, diff)
                checks = checks + msg
        # Now check whether anything failed beyond tolerance
        self.assertFalse(checks, checks)


def main():
    suites = []
    suites.append(unittest.TestLoader().loadTestsFromTestCase(TestVertex))
    suite.append(unittest.TestLoader().loadTestsFromTestCase(TestVertexList))
    if dxf_test:
        suites.append(unittest.TestLoader().loadTestsFromTestCase(TestDXFGeometry))
    alltests = unittest.TestSuite(suites)
    unittest.TextTestRunner(verbosity=2).run(alltests)

# Check whether the script is being excuted by itself
if __name__=='__main__':
    main()
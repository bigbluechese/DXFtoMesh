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

def main():
    suite1 = unittest.TestLoader().loadTestsFromTestCase(TestVertex)
    suite2 = unittest.TestLoader().loadTestsFromTestCase(TestVertexList)
    alltests = unittest.TestSuite([suite1, suite2])
    unittest.TextTestRunner(verbosity=2).run(alltests)

# Check whether the script is being excuted by itself
if __name__=='__main__':
    main()
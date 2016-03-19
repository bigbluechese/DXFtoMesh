#!/usr/bin/python

'''
Code Written by Jeff Peterson (2016)
Takes a DXF file and converts it into segments (which can be arcs or lines) and
vertexes. Currently the code does nothing with circles. 

Requirements:
- Python 2.7x
- dxfgrabber 0.7.5
- matplotlib (if plotting)
- DXFTestSuite.py (if using this for testing)

Usage:
This script can either be used by importing its classes and using them in a
separate script or by running this script directly. The script takes input from
the command line and the usage for which can be found by running this script
with the '-h' flag afterwards specificaly, './DXFtoSegments.py -h'

Typing help(DXFtoSegments) in the interpreter will give you information about
the classes and functions in this module and how they're used.
'''

import dxfgrabber
import argparse
import re
from numbers import Number
if __name__=='__main__':
    import DXFTestSuite
import unittest
import math

def drange(start, stop, n_steps, use_start=True):
    '''
    Generator that creates a decimal range from start to stop using n_steps
    '''
    if n_steps <=1 or type(n_steps) != int:
        raise ValueError('number of steps must be an integer greater than 1')
    if start == stop:
        raise ValueError('Start and stop cannot be the same value')
    incr = (float(stop) - float(start))/(n_steps - 1)
    for i in range(n_steps):
        val = start + i*incr
        if use_start==False and i==0:
            continue
        yield val

def tuple2_check(var):
    '''
    Checkes whether a variable is a tuple of length 2 containing numbers and
    converts numbers to floats if needed

    ARGUMENTS:
    var (tuple)             --  variable to be tested and converted

    RETURNS:
    tup (tuple)             --  tuple of two floats

    RAISES:
    TypeError               --  if vertex is not a tuple
    IndexError              --  if vertex tuple is not length 2
    ValueError              --  if vertex tuple does not contain numbers
    '''
    if type(var) != tuple:
        raise TypeError('variable must be a tuple')
    elif len(var) != 2:
        raise IndexError('tuple must be of length 2')
    elif not(isinstance(var[0], Number)) or not(isinstance(var[1], Number)):
        raise ValueError('tuple must only contain numbers')
    tup = (float(var[0]), float(var[1]))
    return tup

def tuple_string_check(var):
    '''
    Checks whether a given input is a str conversion of a 2-length tuple and
    then ensures that the numbers in the tuple are printed like floats. If a 
    tuple is given instead of a string, it will convert the tuple to the 
    correct string format.

    ARGUMENTS:
    var (str)               --  variable to be tested

    RETURNS:
    str_tup (str)           --  properly formatted string of tuple

    RAISES:
    TypeError               --  if var is not a string or tuple
    ValueError              --  if var is not a length-2 tuple that contained
                                numbers
    '''
    if type(var) == tuple:
        var = str((float(var[0]), float(var[1])))
    elif type(var) != str:
        raise TypeError('variable must be a string conversion of a tuple of two numbers')
    
    # Compile regular expression and check whether tuple is length 2 containing
    # numbers
    tuple_check = re.compile('[(]([-+]?[0-9]*\.?[0-9]+[eE]?[-+]?[0-9]*), ([-+]?[0-9]*\.?[0-9]+[eE]?[-+]?[0-9]*)[)]')
    tuple_match = tuple_check.match(var)
    if tuple_match:
        tuple_float = (float(tuple_match.groups()[0]), float(tuple_match.groups()[1]))
    else:
        raise ValueError('tuple converted to string must have been length 2 and had only numbers')
    
    # Return the new string that is a tuple of two floats
    str_tup = str(tuple_float)
    return str_tup


class Vertex():
    '''
    A class that contains information about the vertex including its coordinates
    and which other vertexes it is connected to.

    ATTRIBUTES:
    x (float)           --  The x-coordinate of the vertex (required)
    y (float)           --  The y-coordinate of the vertex (required)
    id (tuple)          --  A string ID of the vertex that is the tuple pair of
                            coordinates
    connected (list)    --  A list of vertex IDs that are connected to this 
                            vertex

    Notes:
    I thought about having vertexes be movable but decided that this would be
    too complicated to implement with vertex IDs being dependent on the vertex
    positions. If a vertex were to move, it would require changing its name as
    well as changing the connection information in all of the connected
    vertexes. As such, the procedure for "moving" a vertex should be to first
    delete it and then remake it elsewhere with the same connectivity
    information.
    '''
    def __init__(self, coords):
        '''
        ARGUMENTS:
        x (float)           --  The x-coordinate of the vertex (required)
        y (float)           --  The y-coordinate of the vertex (required)

        RAISES:
        TypeError           --  if vertex is not a tuple
        IndexError          --  if vertex tuple is not length 2
        ValueError          --  if vertex tuple does not contain numbers
        '''
        # Make sure input is tuple of length 2 containing numbers and convert
        # numbers to floats
        coords = tuple2_check(coords)

        self.x = coords[0]
        self.y = coords[1]
        self.id = str(coords)
        self.connected = set([])
        self.tuple_check = re.compile('[(]([-+]?[0-9]*\.?[0-9]+[eE]?[-+]?[0-9]*), ([-+]?[0-9]*\.?[0-9]+[eE]?[-+]?[0-9]*)[)]')

    def con(self, vertexID):
        '''
        Records a vertex as connected to the current vertex

        ARGUMENTS:
        vertexID (str)      --  vertexID to connect to the current vertex

        RAISES:
        TypeError           --  if vertexID is not of the form str(tuple)
        RuntimeError        --  if vertexID is the same as this vertex's ID
        '''
        # Check to make sure the ID is a str() of a two-number tuple
        try:
            vertexID = tuple_string_check(vertexID)
        except TypeError:
            raise('vertex to be connected must be a string')
        except ValueError:
            msg = 'Cannot connect vertex \'{}\' because is not the correct form. Vertex IDs must be str() applied to a two-element tuple containing numbers'.format(vertexID)
            raise TypeError(msg)
        if vertexID == self.id:
            raise RuntimeError('Cannot connect a vertex to itself')
        else:
            self.connected.add(vertexID)

    def discon(self, vertexID):
        '''
        Removes a vertexID from the set of connected verticies

        ARGUMENTS:
        vertexID (str)      --  vertexID to disconnect to the current vertex

        RAISES:
        KeyError            --  If vertexID is not part of the set of connected
                                vertecies
        TypeError           --  If vertexID is not of the form str(tuple)
        '''
        # Check to make sure the ID is a str() of a two-number tuple
        try:
            vertexID = tuple_string_check(vertexID)
        except TypeError:
            raise('vertex to be connected must be a string')
        except ValueError:
            msg = 'Cannot disconnect vertex \'{}\' because is not the correct form. Vertex IDs must be str() applied to a two-element tuple containing numbers'.format(vertexID)
            raise TypeError(msg)

        # Raise error if not connected
        try:
            self.connected.remove(vertexID)
        except KeyError:
            raise KeyError('Vertex {} is not connected to {}'.format(vertexID, self.id))
        

class VertexList():
    '''
    A class that provides two ways of keeping track of verticies. The first 
    way is through a python set that contains the tuples that represent vertex
    coordinates. The advantage of a set is that its elements must be unique and
    so you cannot duplicate verticies. The second way is through a dict
    containing the vertex class for each vertex in the geometry. The keys for
    the dict are simple string conversions of the tuple of the vertex
    coordinates. This provides a powerful way to access connectivity information
    about each vertex without looping through lines.

    ATTRIBUTES:
    coordinates (set)       --  A list of all vertex coordinates. Coordinates
                                are given as tuples.
    verticies (dict)        --  A dictionary of all verticies with keys given by
                                the str() of their tuple coodinates. Each item
                                is a vertex class with more specific information
                                about the vertex.
    '''
    def __init__(self):
        self.coordinates = set([]) #set of vertex coordinate tuples
        self.verticies = {} #dictionary of vertex objects

    def add(self, v_coords):
        '''
        Adds a vertex to the list of all verticies

        ARGUMENTS:
        v_coords (tuple)        --  Tuple pair of coordinates for vertex
        '''
        # Make sure input is tuple of length 2 containing numbers and convert
        # numbers to floats
        v_coords = tuple2_check(v_coords)

        # Check current length of vertex set
        initial_len = len(self.coordinates)
        # Try to add the new vertex to the set
        self.coordinates.add(v_coords)
        new_len = len(self.coordinates)

        # Check if the vertex was added. If it wasn't, that means it's already
        # in the set so it shouldn't be added
        if new_len > initial_len:
            self.verticies[str(v_coords)] = Vertex(v_coords)

    def connect(self, vertex1, vertex2):
        '''
        Connect a vertex to another vertex

        ARGUMENTS:
        vertex1 (str)           --  First vertex to be connected
        vertex2 (str)           --  Second vertex to be connected

        RAISES:
        TypeError               --  if vertexes are not tuples converted to
                                    strings.
        KeyError                --  if vertex has not yet been added
        '''
        # Check that input is correct and convert for consistency
        vertex1 = tuple_string_check(vertex1)
        vertex2 = tuple_string_check(vertex2)

        try:
            self.verticies[vertex1].con(vertex2)
        except KeyError as inst:
            raise KeyError('{} is not an existing vertex'.format(vertex1))
        except RuntimeError as inst:
            raise
            
        try:
            self.verticies[vertex2].con(vertex1)
        except KeyError as inst:
            self.verticies[vertex1].discon(vertex2)
            raise KeyError('{} is not an existing vertex'.format(vertex2))

    def disconnect(self, vertex1, vertex2):
        '''
        Disconnect a vertex to another vertex

        ARGUMENTS:
        vertex1 (str)           --  First vertex to be disconnected
        vertex2 (str)           --  Second vertex to be disconnected

        RAISES:
        TypeError               --  if vertexes are not tuples converted to
                                    strings.
        KeyError                --  if vertex has not yet been added
        '''
        # Check that input is correct and convert for consistency
        vertex1 = tuple_string_check(vertex1)
        vertex2 = tuple_string_check(vertex2)

        # Remove vertex 2 from list of connections for vertex 1
        try:
            self.verticies[vertex1].discon(vertex2)
        except KeyError as inst:
            raise
            #raise KeyError('{} is not an existing vertex'.format(vertex1))

        # Remove vertex 1 from list of connections for vertex 2
        try:
            self.verticies[vertex2].discon(vertex1)
        except KeyError:
            raise RuntimeError('Vertex connection not symmetric: {} does not contain connection info for {}'.format(vertex2, vertex1))

    def remove(self, v_coords):
        '''
        Removes a given vertex from the set of verticies

        ARGUMENTS:
        v_coords (tuple)        --  Tuple of coordinates defining the vertex

        RAISES:
        KeyError                -- if v_coords is not an existing vertex
        '''
        # Make sure input is tuple of length 2 containing numbers and convert
        # numbers to floats
        v_coords = tuple2_check(v_coords)

        # Try to remove the vertex and throw exception if it fails
        try:
            self.coordinates.remove(v_coords)
        except KeyError:
            raise KeyError('{} is not an existing vertex'.format(v_coords))

        # Now figure out the connectivity
        connections = self.verticies[str(v_coords)].connected

        # Remove all of the connections between this vertex and others
        for vertexID in connections:
            self.disconnect(str(v_coords), vertexID)

        # Now finally delete the vertex in the dict
        del self.verticies[str(v_coords)]

    def move_vertex(self, v_coords1, v_coords2):
        '''
        Moves a vertex from one position to another by deleting it and then
        recreating it with the same connectivity.

        ARGUMENTS:
        v_coords1 (tuple)       --  Coordinates for the original vertex
        v_coords2 (tuple)       --  New coordinates for the vertex

        RAISES:
        KeyError                --  if v_coords1 or v_coords2 don't exist
        '''
        # Make sure input is tuple of length 2 containing numbers and convert
        # numbers to floats
        v_coords1 = tuple2_check(v_coords1)
        v_coords2 = tuple2_check(v_coords2)

        try:
            connections = self.verticies[str(v_coords1)].connected
        except KeyError:
            raise KeyError('{} is not an existing vertex'.v_coords1)

        # Remove the vertex and its connections
        self.remove(v_coords1)

        # Add the vertex back but at the new coordinates
        self.add(v_coords2)

        # Add back the connections
        for vertexID in connections:
            self.connect(str(v_coords2), vertexID)

class DXFGeometry():
    '''
    Class that first reads a DXF file and then converts the information there
    into a form that is appropriate for making meshes. Specifically, the class
    will take the DXF file, break it into its vertecies and also into the
    respective line segements that exist between vertexes. These line segements
    can either be arcs that are defined in a number of different ways or simple
    straight lines between points.

    ATTRIBUTES:
    dxf (dxfgrabber obj)        --  DXF file read by dxfgrabber module
    verts (VertexList obj)      --  VertexList class containing information
                                    about the verticies in the drawing
    segments (set)              --  Set of segments defined by two points
                                    along with additional information in the
                                    case of arcs/bulges.

    Segments data structure:
    The structure of the segments variable is set of length-2 tuples. Each tuple
    is defined by:
    (1)   a tuple of the endpoint coordinates (each given by a tuple)
    (2)   a list of attributes about the bulge

    The list of bulge attributes is ordered:
    (i)   bulge value
    (ii)  arc beginning angle (relative to horizontal) (radians)
    (iii) arc ending angle (relative to horizontal) (radians)
    (iv)  arc center coordinates
    (v)   radius of arc

    A straight line will contain an empty list for the bulge information
    '''

    def __init__(self, dxf_file):
        '''
        Reads the DXF file and any arguments that may have been passed to the
        class. This could include command-line arguments.
        '''
        self.dxf = dxfgrabber.readfile(dxf_file)
        self.verts = VertexList()
        self.segments = set([])

        # Loop over the entities in the DXF file
        for entity in self.dxf.entities:
            if entity.dxftype == 'LINE':
                self.add_line(entity)
            elif entity.dxftype == 'ARC':
                self.add_arc(entity)
            elif entity.dxftype == 'POLYLINE':
                pass
            elif entity.dxftype == 'CIRCLE':
                pass
            elif entity.dxftype == 'POINT':
                pass

    def add_line(self, dxfentity):
        '''
        Converts a DXF line entity (or entitines) to the proper form and adds
        the information to the list of verticies and to the set of segments

        ARGUMENTS:
        dxfentity (obj)         --  DXF entity or entities from dxfgrabber

        RAISES:

        '''
        try:
            dxftype = dxfentity.dxftype
        except AttributeError:
            pass
        else:
            entities = [dxfentity]

        for entity in entities:
            # Find end-points
            start = (entity.start[0], entity.start[1])
            end = (entity.end[0], entity.end[1])

            # Add verticies
            self.verts.add(start)
            self.verts.add(end)

            # Collect segment data and add to set
            seg = ((start, stop),[])
            self.segments.add(seg)

    def add_arc(self, dxfentity):
        '''
        Converts a DXF arc entity (or entities) to the proper form and adds the 
        information to the list of verticies and to the set of segments. Bulge
        and arc information is also compiled and computed.

        ARGUMENTS:
        entity              --  DXF entity from dxfgrabber

        RAISES:

        '''
        try:
            dxftype = dxfentity.dxftype
        except AttributeError:
            pass
        else:
            entities = [dxfentity]

        for entity in entities:
            # Extract information (again ignoring z-coordinate)
            theta0 = math.radians(entity.startangle)
            theta1 = math.radians(entity.endangle)
            center = (entity.center[0], entity.center[1])
            radius = entity.radius

            # Calculate bulge and start/stop information
            bulge = math.tan((theta1 - theta0)/4)
            start = (radius*math.cos(theta0) + center[0], 
                     radius*math.sin(theta0) + center[1])
            end = (radius*math.cos(theta1) + center[0], 
                   radius*math.sin(theta1) + center[1])

            # Add verticies
            self.verts.add(start)
            self.verts.add(end)

            # Collect segment data and add to set
            seg = ((start, stop), [bulge, theta0, theta1, center, radius])
            self.segments.add(seg)

    def add_polyline(self, dxfentity):
        '''
        '''
        pass

    def rem_reversed(self):
        '''
        '''
        pass



def test_suite(verbose):
    '''Runs the test suite'''
    if verbose:
        verbosity = 2
    else:
        verbosity = 1
    suite1 = unittest.TestLoader().loadTestsFromTestCase(DXFTestSuite.TestVertex)
    suite2 = unittest.TestLoader().loadTestsFromTestCase(DXFTestSuite.TestVertexList)
    alltests = unittest.TestSuite([suite1, suite2])

    unittest.TextTestRunner(verbosity=verbosity).run(alltests)

def main():
    # NOTE: Argument parsing requires Python 2.7x or higher
    # Parses command-line input arguemtns
    help_string = '''Reads a DXF file and then converts it to a suitable form
                     for use in creating computational meshes for CrysMAS and 
                     Cats2D'''
    parser = argparse.ArgumentParser(description=help_string)
    # Specify the DXF file (required)
    help_string = '''Specify the DXF file to convert or specify \'test\' to run
                     the test suite'''
    parser.add_argument('dxf_file', action='store', type=str, 
                        metavar='dxf_file or \'test\'', help=help_string)
    # Create verbose mode
    parser.add_argument('-v', '--verbose', action='store_true')
    # Specify whether information should be output to CrysMAS and optionally specify a new filename
    help_string = '''Turn the DXF file into a CrysMAS .pcs mesh file. A file 
                     name for the .pcs file can optionally be specified'''
    parser.add_argument('-c', '--crysmas', nargs='?', metavar='file',
                        const=True, help=help_string)

    args = parser.parse_args()

    # Specify testing mode from the command line
    if args.dxf_file == 'test':
        testing = True
        test_suite(args.verbose)
    # Otherwise create a DXF geometry object
    else:
        pass

    # Depending on whether verbose mode is enabled, the vprint() function will
    # either print or do nothing. This implementation is faster so that the
    # "if verbose" statement isn't being evaluated each time.
    if args.verbose:
        # Define a verbose print function only when verbose mode is enabled
        def vprint(*args):
            # Basically do the same as print
            for arg in args:
                print arg,
            print
    else:
        # Otherwise, the vprint() function will simply do nothing
        vprint = lambda *a: None




# Check whether the script is being excuted by itself
if __name__=="__main__":
    main()
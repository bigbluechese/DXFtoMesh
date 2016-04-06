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
- DXFTests directory of DXF test files

Usage:
This script can either be used by importing its classes and using them in a
separate script or by running this script directly. The script takes input from
the command line and the usage for which can be found by running this script
with the '-h' flag afterwards specificaly, './DXFtoSegments.py -h'

Typing help(DXFtoSegments) in the interpreter will give you information about
the classes and functions in this module and how they're used.

DEVELOPMENT NOTES:
- In the future it may be advantageous to add segments in batches rather than
individually. This will make it faster to reverse the segments and check them
because you would only have to reverse the segments you're adding. This would
mean compiling a list of segments to be added rather than adding them
individually.
'''

import dxfgrabber
import argparse
import re
from numbers import Number
if __name__=='__main__':
    import DXFTestSuite # Run test suite from this file
import unittest
import math
import numpy as np
from matplotlib import pyplot as plt

def angle360(angle):
    '''Converts an angle from a domain of [-pi, pi] to [0, 2*pi]'''
    if angle < 0:
        new_angle = 2*np.pi + angle
    else:
        new_angle = angle
    return new_angle

def anglespace(start, stop, num=50, endpoint=True, retstep=False):
    '''
    Return evenly spaced angles over a specified angle range in polar
    coordiantes.

    Returns `num` evenly spaced angles calculated over the interval
    [`start`, `stop`] in the counter-clockwise direction. If `start` is greater
    than `stop`, the function will return points that increase from `start`
    until 2*pi and will then restart at zero up until `stop`. A full circle
    would then be specified from 0 to 2*pi.

    ARGUMENTS:
    start : scalar (radians)
        Starting angle for the sequence. Must be within the closed interval 
        [-2*pi, 2*pi].
    stop : scalar (radians)
        End value for the sequence unless `enpoint` is False. In that case, the
        sequence consists of all but the last of num+1 evenly spaced samples
        so that `stop` is excluded. Note that the step size changes when
        `endpoint` is False. Must be within the closed interval 
        [-2*pi, 2*pi].
    num : int [optional]
        Number of samples to generate (default is 50)
    endpoint : bool [optional]
        If True, then `stop` is the last sample. Otherwise it is not included.
        (Default is true)
    retstep : bool [optional]
        If True, return (`samples`, `step`) where `step` is the spacing between
        samples.

    RETURNS:
    samples : ndarray
        `num` equally-spaced samples in the closed interval [start, stop] or
        the half-open interval [start, stop) depending on whether `endpoint` is
        True or False.
    step : float
        Size of spacing between samples
    '''
    # Check if angles are within ranges
    if start < -2*np.pi or start > 2*np.pi:
        raise ValueError('start must be within range [-2*pi, 2*pi]')
    if stop < -2*np.pi or stop > 2*np.pi:
        raise ValueError('stop must be within range [-2*pi, 2*pi]')

    # Make angle always increasing by adding 2*pi to end
    if stop < start:
        end = stop + 2*np.pi
        beg = start
    else:
        beg = start
        end = stop

    # Calculate stepsize
    step = (end - beg)/(num - 1)

    # Remove endpoint if desired
    if not endpoint:
        step = (end - beg)/num
        stop = stop - step

    # Calculate switching index where angle would exceed 360 degrees
    num_switch = int(np.floor(np.around((2*np.pi - start)/step,decimals=12)))+1

    # Create output array from start to zero and then to end or directly to end
    if num_switch < num:
        first_stop = start + (num_switch-1)*step
        first_set = np.linspace(start, first_stop, num_switch)
        second_start = first_stop + step - 2*np.pi
        second_set = np.linspace(second_start, stop, num - num_switch)
        samples = np.concatenate((first_set, second_set))
    else:
        samples = np.linspace(start, stop, num)

    # Return step size if desired
    if retstep:
        return samples, step
    else:
        return samples

def ccw_angle_diff(start, end):
    '''Calculates the difference between two angles in the counter-clockwise
    direction. For example, if the end angle is less than the starting angle,
    this means that the difference will be larger than pi radians'''
    # Check if angles are within ranges
    if start < -2*np.pi or start > 2*np.pi:
        raise ValueError('start must be within range [-2*pi, 2*pi]')
    if end < -2*np.pi or end > 2*np.pi:
        raise ValueError('stop must be within range [-2*pi, 2*pi]')
    # Check whether difference goes through zero
    if end <= start:
        diff = end + 2*np.pi - start
    else:
        diff = end - start
    # Return the difference
    return diff

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
    tol (float)                 --  The tolerance to which all positions are
                                    specified (default 1.0e-08)

    Segments data structure:
    The structure of the segments variable is set of length-2 tuples. Each tuple
    is defined by:
    (1)   a tuple of the endpoint coordinates (each given by a tuple)
    (2)   a tuple of attributes about the bulge

    The tuple of bulge attributes is ordered:
    (i)   bulge value
    (ii)  arc beginning angle (relative to horizontal) (radians)
    (iii) arc ending angle (relative to horizontal) (radians)
    (iv)  arc center coordinates
    (v)   radius of arc

    A straight line will contain an empty tuple for the bulge information
    '''

    def __init__(self, dxf_file, testing=False, verbose=False, tol=1.0e-08):
        '''
        Reads the DXF file and any arguments that may have been passed to the
        class. This could include command-line arguments.

        ARGUMENTS:
        dxf_file (str)          --  Name of DXF file to be loaded

        OPTIONAL ARGUMENTS:
        testing (bool)          --  If testing mode is activated, entities are
                                    NOT added from the DXF file.
                                    (Default = False)
        verbose (bool)          --  When vebose is True, every entity added
                                    will have information about it printed.
                                    (Default = False)
        tol (float)             --  The tolerance to which all positions will be
                                    rounded. This is important when points 
                                    overlap in the DXF file but differ slightly
                                    in their coordinates. Without the tolerance,
                                    these would be treated as separate points
                                    when the user might intend that they are the
                                    same point. This often happens with repeated
                                    patterns.
                                    (Default = 1.0e-08)
        '''
        self.dxf = dxfgrabber.readfile(dxf_file)
        self.dxf_name = dxf_file
        self.verts = VertexList()
        self.segments = set([])
        self.verbose = verbose
        self.tol = tol #Rounds coordinates to this tolerance

        # Add all entities in dxf_file to segments and verts if not testing
        if not testing:
            self.add_entities(self.dxf.entities)

    def vprint(self, *args):
        '''Printing function that is responsive to self.verbose'''
        if self.verbose:
            for arg in args:
                print arg,
            print
    def approx(self, val):
        '''Approximates a value to the given tolerance'''
        decimals = np.abs(np.floor(np.log10(np.abs(self.tol))).astype(int))
        new_val = round(val, decimals)
        if new_val == -0.0:
            new_val = 0.0 #Get rid of negative zero
        #print val, new_val
        return new_val

    def add_dxf(self, dxf_file):
        '''
        Extracts the entities from a dxf file and adds them to the geometry.
        This method is useful for merging two DXF files together.

        ARGUMENTS:
        dxf_file (str)          --  Name of DXF file to be loaded
        '''
        # Consider using a while structure with file to properly handle
        dxf = dxfgrabber.readfile(dxf_file)
        entities = dxf.entities
        dxf = None # Not sure if this is necessary to close file
        self.add_entities(entities)

        return entities


    def add_entities(self, dxfentities):
        '''
        Breaks up DXF entities to extract their segment and vertex information
        and store that information. Also removes reversed duplicates after
        adding the entities.

        ARGUMENTS:
        entities (obj)          --  DXFgrabber entity or entities object to be
                                    added to be stored in segments and verts

        RAISES:
        TypeError               --  if dxftype attribute of the DXF entity is
                                    not 'LINE', 'ARC', 'POLYLINE', 'CIRCLE', or
                                    'POINT'. Currently, only the first three
                                    are handled and the last two are ignored.

        '''
        # Check if object passed is a signle entity or not
        try:
            dxftype = dxfentities.dxftype
        except AttributeError: #Collections don't have dxftype attribute
            entities = dxfentities
        else:
            entities = [dxfentities] #Make into one-element list for looping

        # Loop over entities
        for entity in entities:
            if entity.dxftype == 'LINE':
                self.add_line(entity)
            elif entity.dxftype == 'ARC':
                self.add_arc(entity)
            elif entity.dxftype == 'POLYLINE':
                self.add_polyline(entity)
            elif entity.dxftype == 'CIRCLE':
                pass
            elif entity.dxftype == 'POINT':
                pass
            else:
                raise TypeError('DXF entitiy type, {}, is not supported'.format(entity.dxftype))

    def add_line(self, entity):
        '''
        Converts a DXF line entity (or entitines) to the proper form and adds
        the information to the list of verticies and to the set of segments

        ARGUMENTS:
        entity (obj)            --  DXF entity from dxfgrabber

        RAISES:
        TypeError               --  if entitiy.dxftype is not 'LINE'
        '''
        # Check input
        if entity.dxftype != 'LINE':
            msg = 'dxf entitiy passed was not a LINE but a {}'.format(entity.dxftype)
            raise TypeError(msg)

        # Find end-points
        start = (self.approx(entity.start[0]), self.approx(entity.start[1]))
        end = (self.approx(entity.end[0]), self.approx(entity.end[1]))

        # Add verticies and connect them
        self.verts.add(start)
        self.verts.add(end)
        self.verts.connect(start, end)

        # Collect segment data and add to set
        seg = ((start, end),())
        self.vprint('adding line {}'.format(seg[0]))
        self.segments.add(seg)

    def add_arc(self, entity):
        '''
        Converts a DXF arc entity (or entities) to the proper form and adds the 
        information to the list of verticies and to the set of segments. Bulge
        and arc information is also compiled and computed.

        ARGUMENTS:
        entity (obj)            --  DXF entity from dxfgrabber

        RAISES:
        TypeError               --  if entitiy.dxftype is not 'ARC'
        '''
        # Check input
        if entity.dxftype != 'ARC':
            msg = 'dxf entitiy passed was not an ARC but a {}'.format(entity.dxftype)
            raise TypeError(msg)

        # Extract information (again ignoring z-coordinate)
        start_angle = np.radians(entity.startangle)
        end_angle = np.radians(entity.endangle)
        center = (entity.center[0], entity.center[1])
        radius = entity.radius

        # Calculate bulge and start/stop information
        theta = ccw_angle_diff(start_angle, end_angle)
        bulge = np.tan(theta/4)
        start = (self.approx(radius*np.cos(start_angle) + center[0]), 
                 self.approx(radius*np.sin(start_angle) + center[1]))
        end = (self.approx(radius*np.cos(end_angle) + center[0]), 
               self.approx(radius*np.sin(end_angle) + center[1]))

        # Add verticies and connect them
        self.verts.add(start)
        self.verts.add(end)
        self.verts.connect(start, end)

        # Collect segment data and add to set
        seg = ((start, end), (bulge, start_angle, end_angle, center, radius))
        self.vprint('adding arc {}'.format(seg[0]))
        self.segments.add(seg)

    def add_polyline(self, entity):
        '''
        Converts a DXF polyline entity into segments while transforming bulges
        into arcs (defined by a center of curvature and radius of curvature) and
        storing the vertex and segment information.

        ARGUMENTS:
        entity (obj)            --  DXF entity from dxfgrabber
        
        RAISES:
        TypeError               --  if entitiy.dxftype is not 'POLYLINE'
        '''
        # Check input
        if entity.dxftype != 'POLYLINE':
            msg = 'dxf entitiy passed was not a POLYLINE but a {}'.format(entity.dxftype)
            raise TypeError(msg)

        self.vprint('Breaking up this polyline into segments:\n{}'.format(entity))
        # Loop through the points in the polyline
        for i, point in enumerate(entity.points):
            # Add the current point
            start = (self.approx(point[0]), self.approx(point[1]))
            self.verts.add(start)
            try:
                # Add the next point if it exists
                next_point = entity.points[i+1]
            except IndexError:
                # Next point DOESN'T exist therefore this is the end of the
                # polyline
                if entity.is_closed:
                    # If polyline is closed, connect the last point (the current
                    # point) back to the first point
                    first_point = entity.points[0]
                    end = (self.approx(first_point[0]), 
                           self.approx(first_point[1]))
                    self.verts.connect(start, end)
                else:
                    # Otherwise the polyline is open so all segments have been
                    # added already
                    self.vprint('\tThis polyline is not closed!')
                    break
            else:
                # The next point DOES exist so add it and connect it to the
                # current point
                end = (self.approx(next_point[0]), self.approx(next_point[1]))
                self.verts.add(end)
                # Connect the two points
                #print start, end
                self.verts.connect(start, end)

            # Check whether there is a bulge in this segment
            if  entity.bulge[i] != 0:
                # Convert bulge information to arc and store information
                bulge = entity.bulge[i]
                # Distance between points
                d = np.sqrt((start[0] - end[0])**2 + (start[1] - end[1])**2)
                # Angle between points from center
                theta = 4*np.arctan(bulge)
                # Radius of circle making arc
                radius = d/2/np.sin(abs(theta)/2)
                # Find angle of segment relative to x axis
                alpha = np.arctan2(end[1]-start[1], end[0]-start[0])
                # Find angle between segment and radius (obtuse is negative)
                beta = (np.pi/2 - abs(theta)/2)*\
                           (np.pi - theta)/abs(np.pi - theta)
                # Angle to radius vector from x-axis is then the sum of alpha
                # and beta
                gamma = alpha - beta
                # Gamma angle and radius describe the vector pointing from the
                # start point to the center
                center = (radius*np.cos(gamma)+start[0],
                          radius*np.sin(gamma)+start[1])
                # Now compute start and stop angles relative to horizontal in
                # a counter-clockwise sense
                start_angle = angle360(np.arctan2(start[1]-center[1], 
                                                  start[0]-center[0]))
                end_angle = angle360(np.arctan2(end[1]-center[1],
                                                end[0]-center[0]))

                # Compile all bulge/arc information and add it to segments
                seg = ((start, end), (bulge, start_angle, end_angle, center,
                                      radius))
                self.vprint('\tadding arc {}'.format(seg[0]))
                self.segments.add(seg)

            # Segment is a straight line
            else:
                seg = ((start, end), ())
                # Add info to segments
                self.vprint('\tadding line {}'.format(seg[0]))
                self.segments.add(seg)

    def rem_reversed(self):
        '''
        Looks at the current set of segments and removes repeat segments that
        exist as reversals of currently existing segments. This is accomplished
        by reversing every segment (inclduing recalculating arc/bulge info) and
        then trying to remove that segment from the segment set.
        '''

        pruned_segs = self.segments.copy()
        for seg in self.segments.copy():
            # First reverse the segment
            rev_seg_coords = (seg[0][1], seg[0][0])
            if seg[1] == ():
                rev_seg_info = ()
            else:
                # Calculate reversed arc/bulge information
                bulge, start_angle, end_angle, center, radius = seg[1]
                rev_bulge = -bulge
                rev_s_angle = end_angle
                rev_e_angle = start_angle
                rev_center = center
                rev_radius = radius
                rev_seg_info = (rev_bulge, rev_s_angle, rev_e_angle, rev_center,
                    rev_radius)
            # Create the reversed segment
            rev_seg = (rev_seg_coords, rev_seg_info)
            # Now determine if this segment already exists in the set
            try:
                pruned_segs.remove(rev_seg)
            except KeyError:
                continue #It's not in the segment list so move on
        self.segments = pruned_segs
        return pruned_segs

    def display(self):
        '''
        Plots the current geometry but does not display it. In order to display
        the plot, a plt.show() call must be made after the display() method is
        used with plt of course being the matplotlib.pyplot module.
        '''
        def arc_points(bulge, center, radius, startangle, endangle):
            '''
            Creates a series of points that form an arc. Outputs series of x
            points and y points that form arc.
            '''
            if bulge > 0:
                angles = np.linspace(startangle, endangle, num=50)
            else:
                angles = np.linspace(endangle, startangle, num=50)
            x_vals = radius*np.cos(angles) + center[0]
            y_vals = radius*np.sin(angles) + center[1]
            return (x_vals, y_vals)

        def arc_points2(bulge, center, radius, startangle, endangle):
            '''
            Creates a series of points that form an arc. Outputs series of x
            points and y points that form arc. Correctly accounts for the fact
            that the start angle could be larger than the end angle.
            '''
            # If the bulge is less than 0, arc must be defined clockwise
            if bulge > 0:
                angles = anglespace(startangle, endangle, num=50)
            else:
                angles = anglespace(endangle, startangle, num=50)
            x_vals = radius*np.cos(angles) + center[0]
            y_vals = radius*np.sin(angles) + center[1]
            return (x_vals, y_vals)


        def seg_plot(segment):
            '''
            Creates information for a plot from a segment
            '''
            # If segment is a straight line, this will evaluate true
            if not segment[1]:
                x_vals, y_vals = zip(*segment[0])
            # Otherwise extract the arc information and generate x,y points
            else:
                bulge = segment[1][0]
                startangle = segment[1][1]
                endangle = segment[1][2]
                center = segment[1][3]
                radius = segment[1][4]
                x_vals, y_vals = arc_points2(bulge, center, radius, startangle,
                                            endangle)
            return (x_vals, y_vals)

        # Set up the plot space
        fig1, ax1 = plt.subplots()

        # Plot the lines and arcs
        for seg in self.segments:
            x, y = seg_plot(seg)
            ax1.plot(x,y, 'b')

        # Plot vertex locations
        x_coords, y_coords = zip(*self.verts.coordinates)
        ax1.plot(x_coords, y_coords, 'ks')

        # Create the plot
        plt.axis('scaled')
        plt.draw() # Plot must be shown to be visible so after calling the

def test_suite(verbose=False, dxf_test=False):
    '''Runs the test suite'''
    if verbose:
        verbosity = 2
    else:
        verbosity = 1
    suites = []
    suites.append(unittest.TestLoader().loadTestsFromTestCase(DXFTestSuite.TestVertex))
    suites.append(unittest.TestLoader().loadTestsFromTestCase(DXFTestSuite.TestVertexList))
    if dxf_test:
        suites.append(unittest.TestLoader().loadTestsFromTestCase(DXFTestSuite.TestDXFGeometry))
        dxf1 = DXFGeometry('./DXFTests/DXFTest1.dxf')
        dxf1.display()
        dxf2 = DXFGeometry('./DXFTests/DXFTest2.dxf')
        dxf2.display()
        dxf3 = DXFGeometry('./DXFTests/DXFTest3.dxf')
        dxf3.display()
        dxf4 = DXFGeometry('./DXFTests/DXFTest4.dxf')
        dxf4.display()
    alltests = unittest.TestSuite(suites)
    unittest.TextTestRunner(verbosity=verbosity).run(alltests)
    # Show the plotted DXF file if available
    plt.show()

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
    # Skip DXF tests if option is passed
    help_string = '''Skips the DXF tests if testing mode is activated'''
    parser.add_argument('--nodxf', action='store_true', help=help_string)

    args = parser.parse_args()

    # Specify testing mode from the command line
    if args.dxf_file == 'test':
        testing = True
        test_suite(verbose=args.verbose, dxf_test=not(args.nodxf))
    # Otherwise create a DXF geometry object
    else:
        dxf = DXFGeometry(args.dxf_file, verbose=args.verbose)


# Check whether the script is being excuted by itself
if __name__=="__main__":
    main()
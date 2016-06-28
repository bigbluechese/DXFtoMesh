#!/usr/bin/env python
'''
Code Written by Jeff Peterson (2016)

The DXFGeometry object in this file is the primary class. It can be used to read
a DXF file and then convert the information into a form that can then be used
in the construction of a computational mesh. This module should be imported into
another script and the classes used there.

DEVELOPMENT NOTES:
- In the future it may be advantageous to add segments in batches rather than
individually. This will make it faster to reverse the segments and check them
because you would only have to reverse the segments you're adding. This would
mean compiling a list of segments to be added rather than adding them
individually.
- Segments maybe should be made into a separate class so that reversing the
segment list and looking for duplicates can be completely avoided. In this way,
a reversed segment and a forward segment will be the same and thus wouldn't
be duplicated in the set.
'''

from __future__ import print_function
import dxfgrabber
import argparse
import re
from numbers import Number
import math
import numpy as np
import os
import pickle
import pkg_resources
from matplotlib import pyplot as plt
from HelperFunctions import angle360, anglespace, approx, bulge_to_arc, \
                            ccw_angle_diff, tuple2_check
import lineintersect

class Vertex(object):
    '''
    A class that contains information about the vertex including its coordinates
    and which other vertexes it is connected to.

    ATTRIBUTES:
    x (float)           --  The x-coordinate of the vertex (required)
    y (float)           --  The y-coordinate of the vertex (required)
    id (tuple)          --  A string ID of the vertex that is the tuple pair of
                            coordinates
    connections (list)  --  A list of vertex IDs that are connected to this 
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
        coords (tuple)      --  The x- and y-coordinates of the vertex

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
        self.id = coords
        self.connections = set([])

    def __repr__(self):
        return 'Vertex({}, {})'.format(self.x, self.y)

    def con(self, vertexID):
        '''
        Records a vertex as connected to the current vertex

        ARGUMENTS:
        vertexID (str)      --  vertexID to connect to the current vertex

        RAISES:
        TypeError           --  if vertexID is not a tuple
        RuntimeError        --  if vertexID is the same as this vertex's ID
        '''
        # Check to make sure the ID is a two-number tuple
        try:
            vertexID = tuple2_check(vertexID)
        except (TypeError, IndexError) as e:
            raise('''vertex to be connected must be a length 2 tuple not 
                  {}'''.format(vertexID))
        except ValueError:
            msg = '''Cannot connect vertex \'{}\' because is does not contain 
            numbers'''.format(vertexID)
            raise TypeError(msg)
        if vertexID == self.id:
            raise RuntimeError('Cannot connect a vertex to itself')
        else:
            self.connections.add(vertexID)

    def discon(self, vertexID):
        '''
        Removes a vertexID from the set of connected vertices

        ARGUMENTS:
        vertexID (str)      --  vertexID to disconnect to the current vertex

        RAISES:
        KeyError            --  If vertexID is not part of the set of connected
                                vertices
        TypeError           --  If vertexID is not a tuple
        '''
        # Check to make sure the ID is a two-number tuple
        try:
            vertexID = tuple2_check(vertexID)
        except (TypeError, IndexError) as e:
            raise('''vertex to be connected must be a length 2 tuple not 
                  {}'''.format(vertexID))
        except ValueError:
            msg = '''Cannot connect vertex \'{}\' because is does not contain 
            numbers'''.format(vertexID)
            raise TypeError(msg)

        # Raise error if not connected
        try:
            self.connections.remove(vertexID)
        except KeyError:
            raise KeyError('Vertex {} is not connected to {}'.format(vertexID, self.id))
        

class VertexList(object):
    '''
    A class that provides two ways of keeping track of vertices. The first 
    way is through a python set that contains the tuples that represent vertex
    coordinates. The advantage of a set is that its elements must be unique and
    so you cannot duplicate vertices. The second way is through a dict
    containing the vertex class for each vertex in the geometry. The keys for
    the dict are simple string conversions of the tuple of the vertex
    coordinates. This provides a powerful way to access connectivity information
    about each vertex without looping through lines.

    DEVELOPER NOTES:
    Is the coordinates set even relevant? Consider checking the dictionary prior
    to any addition and using that instead. This would also allow the use of an
    inheritence so that the VertexList could behave as a VertexDict.

    ATTRIBUTES:
    coordinates (set)       --  A list of all vertex coordinates. Coordinates
                                are given as tuples.
    vertices (dict)        --  A dictionary of all vertices with keys given by
                                their tuple coodinates. Each item is a vertex
                                class with more specific information about the
                                vertex.
    '''
    def __init__(self):
        self.coordinates = set([]) #set of vertex coordinate tuples
        self.vertices = {} #dictionary of vertex objects

    def add(self, v_coords):
        '''
        Adds a vertex to the list of all vertices

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
            self.vertices[v_coords] = Vertex(v_coords)

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
        vertex1 = tuple2_check(vertex1)
        vertex2 = tuple2_check(vertex2)

        try:
            self.vertices[vertex1].con(vertex2)
        except KeyError as inst:
            raise KeyError('{} is not an existing vertex'.format(vertex1))
        except RuntimeError as inst:
            raise
            
        try:
            self.vertices[vertex2].con(vertex1)
        except KeyError as inst:
            self.vertices[vertex1].discon(vertex2)
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
        vertex1 = tuple2_check(vertex1)
        vertex2 = tuple2_check(vertex2)

        # Remove vertex 2 from list of connections for vertex 1
        try:
            self.vertices[vertex1].discon(vertex2)
        except KeyError as inst:
            raise

        # Remove vertex 1 from list of connections for vertex 2
        try:
            self.vertices[vertex2].discon(vertex1)
        except KeyError:
            raise RuntimeError('''Vertex connection not symmetric: {} does not
                    contain connection info for {}'''.format(vertex2, vertex1))

    def remove(self, v_coords):
        '''
        Removes a given vertex from the set of vertices

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
        connections = self.vertices[v_coords].connections.copy()

        # Remove all of the connections between this vertex and others
        for vertexID in connections:
            self.disconnect(v_coords, vertexID)

        # Now finally delete the vertex in the dict
        del self.vertices[v_coords]

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
            connections = self.vertices[v_coords1].connections
        except KeyError:
            raise KeyError('{} is not an existing vertex'.v_coords1)

        # Remove the vertex and its connections
        self.remove(v_coords1)

        # Add the vertex back but at the new coordinates
        self.add(v_coords2)

        # Add back the connections
        for vertexID in connections:
            self.connect(v_coords2, vertexID)

class DXFGeometry(object):
    '''
    DEVELOPER NOTES:
        - Need to add a method to move a vertex

    Class that first reads a DXF file and then converts the information there
    into a form that is appropriate for making meshes. Specifically, the class
    will take the DXF file, break it into its vertices and also into the
    respective line segements that exist between vertexes. These line segements
    can either be arcs that are defined in a number of different ways or simple
    straight lines between points.

    ATTRIBUTES:
    dxf (dxfgrabber obj)        --  DXF file read by dxfgrabber module
    verts (VertexList obj)      --  VertexList class containing information
                                    about the vertices in the drawing
    segments (set)              --  Set of segments defined by two points
                                    along with additional information in the
                                    case of arcs/bulges.
    tol (float)                 --  The tolerance to which all positions are
                                    specified (default 1.0e-08)
    work_dir (str)              --  Working directory for the CrysMAS file

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
        # Make sure the right dxfgrabber is installed
        self.dxfgrabber_version_check()

        self.dxf_path = dxf_file
        self.verts = VertexList()
        self.segments = set([])
        self.tol = tol #Rounds coordinates to this tolerance
        self.testing = testing
        no_file = False

        # Turn on/off verbose printing
        self.__verbose__ = verbose
        if self.__verbose__:
            self.turn_on_verbose()
        else:
            self.turn_off_verbose()

        # Turn off verbose mode while testing
        if self.testing:
            self.turn_off_verbose()

        # Extract path information
        self.work_dir, self.dxf_name = os.path.split(os.path.splitext(self.dxf_path)[0])

        # If testing, don't worry about file name
        try:
            self.dxf = dxfgrabber.readfile(dxf_file)
        except IOError:
            if not self.testing:
                raise
            else:
                no_file = True

        # Add all entities in dxf_file to segments and verts if not testing
        if not no_file:
            self.add_entities(self.dxf.entities)
            self.rem_reversed()

    def __repr__(self):
        return 'DXFGeometry object created from {}.dxf'.format(self.dxf_name)

    def dxfgrabber_version_check(self):
        '''
        Checks if the current DXFgrabber version is acceptable. Currently
        version 0.8.0 is incompatible.
        '''
        version = pkg_resources.get_distribution('dxfgrabber').version
        if version == '0.8.0':
            raise ImportError('dxfgrabber version should not be 0.8.0')
        elif cmp('0.7.5', version) > 0:
            print('''WARNING: dxfgrabber versions later than 0.8.0 have not been
                     tested... use with caution''')

    def turn_on_verbose(self):
        '''Turns on verbose printing'''
        self.vprint = lambda *args: print(*args) # Standard print function

    def turn_off_verbose(self):
        '''Turns off verbose printing'''
        self.vprint = lambda *args: None # Return nothing

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
            elif entity.dxftype == 'POLYLINE' or entity.dxftype == 'LWPOLYLINE':
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
        the information to the list of vertices and to the set of segments

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
        start = (approx(entity.start[0], tol=self.tol), approx(entity.start[1], tol=self.tol))
        end = (approx(entity.end[0], tol=self.tol), approx(entity.end[1], tol=self.tol))

        # Add vertices and connect them
        self.verts.add(start)
        self.verts.add(end)
        self.verts.connect(start, end)

        # Collect segment data and add to set
        initial_len = len(self.segments)
        seg = ((start, end),())
        self.vprint('adding line {}'.format(seg[0]))
        self.segments.add(seg)
        if len(self.segments) == initial_len:
            self.vprint('\tSegment already exists... skipped')

    def add_arc(self, entity):
        '''
        Converts a DXF arc entity (or entities) to the proper form and adds the 
        information to the list of vertices and to the set of segments. Bulge
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
        center = (approx(entity.center[0], tol=self.tol), approx(entity.center[1], tol=self.tol))
        radius = approx(entity.radius, tol=1.0e-12)

        # Calculate bulge and start/stop information
        theta = ccw_angle_diff(start_angle, end_angle)
        bulge = np.tan(theta/4)
        start = (approx(radius*np.cos(start_angle) + center[0], tol=self.tol), 
                 approx(radius*np.sin(start_angle) + center[1], tol=self.tol))
        end = (approx(radius*np.cos(end_angle) + center[0], tol=self.tol), 
               approx(radius*np.sin(end_angle) + center[1], tol=self.tol))

        # Add vertices and connect them
        self.verts.add(start)
        self.verts.add(end)
        self.verts.connect(start, end)

        # Collect segment data and add to set
        initial_len = len(self.segments)
        seg = ((start, end), (bulge, start_angle, end_angle, center, radius))
        self.vprint('adding arc {}'.format(seg[0]))
        self.segments.add(seg)
        if len(self.segments) == initial_len:
            self.vprint('\tSegment already exists... skipped')

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
        if entity.dxftype != 'POLYLINE' and entity.dxftype != 'LWPOLYLINE':
            msg = 'dxf entitiy passed was not a POLYLINE but a {}'.format(entity.dxftype)
            raise TypeError(msg)

        self.vprint('Breaking up this polyline into segments:\n{}'.format(entity))
        # Loop through the points in the polyline
        for i, point in enumerate(entity.points):
            # Add the current point
            start = (approx(point[0], tol=self.tol), 
                     approx(point[1], tol=self.tol))
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
                    end = (approx(first_point[0], tol=self.tol), 
                           approx(first_point[1], tol=self.tol))
                    self.verts.connect(start, end)
                else:
                    # Otherwise the polyline is open so all segments have been
                    # added already
                    self.vprint('\tThis polyline is not closed!')
                    break
            else:
                # The next point DOES exist so add it and connect it to the
                # current point
                end = (approx(next_point[0], tol=self.tol), 
                       approx(next_point[1], tol=self.tol))
                self.verts.add(end)
                # Connect the two points
                #print start, end
                self.verts.connect(start, end)

            # Find number of segments before adding this one
            initial_len = len(self.segments)

            # Check whether there is a bulge in this segment
            if  entity.bulge[i] != 0:
                # Convert bulge information to arc and store information
                bulge = entity.bulge[i]
                # Distance between points
                d = np.sqrt((start[0] - end[0])**2 + (start[1] - end[1])**2)
                # Angle between points from center
                theta = 4*np.arctan(bulge)
                # Radius of circle making arc
                radius = approx(d/2/np.sin(abs(theta)/2), tol=1e-12)
                # Find angle of segment relative to x axis
                alpha = np.arctan2(end[1]-start[1], end[0]-start[0])
                # beta = (np.pi/2 - abs(theta)/2)*(np.pi - abs(theta))/abs(np.pi - abs(theta))
                # # Angle to radius vector from x-axis is then the sum of alpha
                # # and beta
                # gamma = alpha + beta
                if bulge > 0:
                    # Find angle between segment and radius. Beta is negative if
                    # theta is greater than pi
                    beta = np.pi/2 - theta/2
                    # Angle to radius vector from x-axis is then the SUM of
                    # alpha and beta
                    gamma = alpha + beta
                else:
                    # Find angle between segment and radius. Beta is negative if
                    # theta is greater than pi
                    beta = np.pi/2 + theta/2
                    # Angle to radius vector from x-axis is then the DIFFERENCE
                    # between alpha and beta
                    gamma = alpha - beta
                # Gamma angle and radius describe the vector pointing from the
                # start point to the center
                center = (approx(radius*np.cos(gamma)+start[0], tol=self.tol),
                          approx(radius*np.sin(gamma)+start[1], tol=self.tol))
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
                if len(self.segments) == initial_len:
                    self.vprint('\tSegment already exists... skipped')

            # Segment is a straight line
            else:
                seg = ((start, end), ())
                # Add info to segments
                self.vprint('\tadding line {}'.format(seg[0]))
                self.segments.add(seg)
                if len(self.segments) == initial_len:
                    self.vprint('\tSegment already exists... skipped')

    def reverse_seg(self, seg):
        '''
        Reverses a given segment including the bulge information
        '''
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
        return rev_seg

    def rem_reversed(self):
        '''
        Looks at the current set of segments and removes repeat segments that
        exist as reversals of currently existing segments. This is accomplished
        by reversing every segment and then trying to remove the original, non-
        reversed segment if the reversed segment also exists.
        '''
        # Convert segment set to list
        pruned_segs = self.segments.copy()
        for seg in self.segments:
            rev_seg = self.reverse_seg(seg)
            # Check if the reversed segment is in the list of segments
            if rev_seg in pruned_segs:
                # Remove the non-reversed segment
                pruned_segs.remove(seg)
            else:
                continue
        self.segments = pruned_segs
        if not self.testing:
            self.vprint('Reversed segments have been removed from {}'.format(self.dxf_path))
        return pruned_segs

    def move_vertex(self, old_coords, new_coords):
        '''
        Moves a vertex belonging to a line. At the moment, only straight lines
        are supported.

        ARGUMENTS:
        old_coords (tuple)  --  Coordinates of the old vertex expressed as a
                                length-2 tuple
        new_coords (tuple)  --  Coordiantes of the old vertex expressed as arc
                                length-2 tuple

        RAISES:
        RuntimeError        --  if a line segment either appears to not exist or
                                if the segment contains bulge information

        DEVELOPER NOTES:
        If support for bulges were to be added, one way would be to add a
        segments dict where the vertices were used as a key and the information
        was a tuple of bulge information entries. More than one entry could
        exist for a given set of vertices since one can have multiple lines
        between two points with different bulge values. This would allow bulge
        information to be modified by moving the vertex although it would be
        tricky to rationally choose how to modify the bulge information.
        '''
        # Make sure the old_coords are actually a vertex
        try:
            old_vert = self.verts.vertices[old_coords]
        except KeyError:
            print('Vertex cannot be moved because the vertex does not exist!')
            raise

        # Now find segment information that corresponds to this vertex
        old_lines = []
        new_lines = []
        for vert in old_vert.connections:
            if ((old_coords, vert), ()) in self.segments:
                old_lines.append(((old_coords, vert), ()))
                new_lines.append(((new_coords, vert), ()))
            elif ((vert, old_coords), ()) in self.segments:
                old_lines.append(((vert, old_coords), ()))
                new_lines.append(((vert, new_coords), ()))
            else:
                raise RuntimeError('''The segment ({}, {}) contains bulge
                        information. Segments with bulges cannot be moved at this
                        time.'''.format(old_coords, vert))

        # Move the vertex in the vertex list
        self.verts.move_vertex(old_coords, new_coords)
        # Recreate the line segments
        for old, new in zip(old_lines, new_lines):
            self.segments.remove(old)
            self.segments.add(new)

    def find_intersections(self):
        '''
        Finds the intersections between line segments. These intersections can
        be of four forms:
        1) two lines simply intersect and have different slopes
        2) two lines share a common end-point
        3) one endpoint lies on the other line
        4) both lines are colinear and overlap

        or no intersection occurs. A sweep-line algorithm is used to first order
        the vertices and then sweep through the verticies in order from left to
        right. If a right-endpoint of a line causes that line to switch
        positions in the currently active lines (ordered from bottom to top),
        an intersection has likely occured and must be tested for.
        '''
        pass

    def cats2d_convert(self, invert_coords=True, len_scale=None):
        '''
        Converts the data to a form that can be used in creating a Cats2D mesh.

        OPTIONAL ARGUMENTS:
        invert_coords (bool)--  When evalatues to True, the coordinates for
                                every point will be swapped. This facilitates
                                the fact that Cats2D operates in a geometry
                                where y maps to r and x maps to z when a CAD
                                drawing will likely map the opposite.
                                (Default = True)
        len_scale (float)   --  If specified, all coordinates in the DXF file
                                will be divided by this value. It is up to the
                                user to ensure that the correct units are used.

        RETURNS:
        v_coords (list)     --  List of vertex coordinates where each coordinate
                                is given by a tuple.
        edges (list)        --  A list of tuples where each tuple contains two
                                indicies identifying which vertices in v_coords
                                make up the edge
        bulges (list)       --  If any edges have bulges, the bulge information
                                is saved and indexed to the edges in the form
                                (edge_index, (bulge_info)).

        WARNINGS:
        UserWarning         --  if bulge information is passed to this function.
                                Currently, the bulge information is not output
                                to Cats2D and may or may not be destroyed.
        '''
        # Find the intersections

        # Do something for each type of intersection

        # Create the v_coords list
        v_coords = list(self.verts.coordinates)
        # Now create a dictionary to associate coordinates with an index
        v_dict = dict(zip(v_coords, range(len(v_coords))))

        # Convert segments to edges by looping and converting vertex coordinates
        # to vertex indicies
        edges = []
        bulges = []
        for seg in self.segments:
            start_vert = v_dict[seg[0][0]]
            end_vert = v_dict[seg[0][1]]
            edges.append((start_vert, end_vert))
            # If the segment is a bulge, save the information in bulges
            if seg[1]:
                i = len(edges) - 1 #Find edge index
                bulges.append((i, seg[1]))
                print('''WARNING: Segment {} contains a bulge.\n\tBulge Info: {}
                      '''.format(seg[0], seg[1]))

        # Swap the x and y coordinates by default
        if invert_coords:
            x_coords, y_coords = zip(*v_coords)
            v_coords = zip(y_coords, x_coords)
            self.vprint('Swapping x and y coordinates for export to Cats2D')

        # Non-dimensionalize the coordiantes if needed
        if len_scale:
            len_scale = float(len_scale)
            v_coords = [(x/len_scale, y/len_scale) for x, y in v_coords]

        # Now return information
        return v_coords, edges, bulges

    def output_to_crysmas(self, dxf_units='mm', f_name=None):
        '''
        Outputs the current geometry to a CrysMAS mesh file. Will convert from
        the dxf_units into meters by the appropriate conversion factor. Bulges
        are converted into straight lines.

        NOTE: There is a chance that CrysMAS will have issues with Windows
        newline characters. If this is the case, run this code on a Unix
        platform to create the .pcs file.

        OPTIONAL ARGUMENTS:
        dxf_units (str)     --  Specify the units of the DXF file as either
                                'mm', 'cm', 'm', 'in' or 'ft'. The code will
                                then convert the dimensions into meters for
                                CrysMAS (Default='mm').
        f_name (str)        --  Output .pcs file name. If none is specified,
                                the name of the DXF file is used for the .pcs
                                file as well.

        RAISES:
        TypeError           --  if units are given that aren't any of the above
        '''
        # Create scaling factor for units (converting to m)
        units = {'m':1, 'cm':0.01, 'mm':0.001, 'in':2.54/100., 'ft':12*2.54/100}
        try:
            scale_factor = units[dxf_units]
        except IndexError:
            msg = 'dxf units must be specified in {}'.format(', '.join(units.keys()))
            raise TypeError(msg)

        # Open the CrysMAS .pcs file
        if f_name == None:
            f_name = self.dxf_name
        crys_path = os.path.join(self.work_dir, f_name+'.pcs')
        crys_file = open(crys_path, 'w')

        # Find extra information for CrysMAS
        x_vals, y_vals = zip(*self.verts.coordinates)
        x_max = max(x_vals)*scale_factor
        x_min = min(x_vals)*scale_factor
        y_max = max(y_vals)*scale_factor
        y_min = min(y_vals)*scale_factor
        num_points = len(self.verts.coordinates)
        num_lines = len(self.segments)

        # Write the max/min header
        line = '{} {} {} {}\n'.format(x_min, x_max, y_min, y_max)
        crys_file.write(line)

        # Next write the number of points
        line = '{} points\n'.format(num_points)
        crys_file.write(line)

        # Create dictionary for matching vertex coordinates to indicies
        v_dict = {}

        # Loop through vertices, assign indicies, and write to file
        for i, v in enumerate(self.verts.coordinates):
            v_scaled = (v[0]*scale_factor, v[1]*scale_factor)
            v_dict[v_scaled] = i+1 #Index is CrysMAS (i.e. starts at 1)
            line = '{} {} {}\n'.format(i+1, v_scaled[0], v_scaled[1])
            crys_file.write(line)

        # Write the number of lines
        line = '{} lines\n'.format(num_lines)
        crys_file.write(line)

        # Loop through segments, assign vertices, and look up vertices
        for i, seg in enumerate(self.segments):
            start_coords = (seg[0][0][0]*scale_factor, seg[0][0][1]*scale_factor)
            end_coords = (seg[0][1][0]*scale_factor, seg[0][1][1]*scale_factor)
            start_vert = v_dict[start_coords]
            end_vert = v_dict[end_coords]
            line = '{} {} {} 0\n'.format(i, start_vert, end_vert)
            crys_file.write(line)

        # Close the file
        crys_file.close()

        print('Saving to CrysMAS geometry {}...'.format(crys_path))


    def display(self):
        '''
        Plots the current geometry but does not display it. In order to display
        the plot, a plt.show() call must be made after the display() method is
        used with plt of course being the matplotlib.pyplot module.
        '''
        def arc_points(bulge, center, radius, startangle, endangle):
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
                x_vals, y_vals = arc_points(bulge, center, radius, startangle,
                                            endangle)
            return (x_vals, y_vals)

        # Set up the plot space
        self.fig, self.ax = plt.subplots()

        # Plot the lines and arcs
        for seg in self.segments:
            x, y = seg_plot(seg)
            self.ax.plot(x,y, 'b')

        # Plot vertex locations
        x_coords, y_coords = zip(*self.verts.coordinates)
        self.ax.plot(x_coords, y_coords, 'ks')

        # Create the plot
        plt.axis('scaled')
        plt.draw() # Plot must be shown to be visible so after calling the

        print('Geometry for {} is queued for display...'.format(self.dxf_path))

    def pickle(self, f_name=None):
        '''
        Outputs the current geometry to a user-specified pickle file
        '''
        if f_name == None:
            f_name = self.dxf_name+'.pickle'

        # Dump the pickle to file
        fo = open(os.path.join(self.work_dir, f_name), 'wb')
        pickle.dump(self, fo)
        fo.close()

def main():
    print('''This file contains classes that are used to create a DXFGeometry
    object from a DXF file for then creating a computational mesh in either
    CrysMAS or Cats2D. Please run the MeshMaker.py file for usage
    information''')

# Check whether the script is being excuted by itself
if __name__=="__main__":
    main()
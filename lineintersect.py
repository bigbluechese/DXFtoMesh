#!/usr/bin/env python

'''
Calculates the intersections points of a bunch of lines and then breaks those
lines up into segments such that no lines intersect and instead only share
end points
'''

import math
from sortedcontainers import SortedDict

class active_vertex(object):
    '''
    Active vertex class for use in finding intersections between lines. An
    active vertex is one that the sweep line has passed over and that has
    connected vertices over which the sweep line has NOT yet passed.

    ATTRIBUTES:
    x (float)           --  x-coordinate of vertex
    y (float)           --  y-coordiante of vertex
    connections (list)  --  list of tuples that represented connected vertices
    rem_req (int)       --  number of remove requests that have been issued
                            against the vertex

    RAISES:
    RuntimeError        --  if the vertex has no connections and thus is useless
                            for defining a line
    '''
    def __init__(self, vert_obj):
        '''
        ARGUMENTS:
        vert_obj (vertex object)
        --  Vertex object containing information about the vertex with the
            following attributes:
                - x (float)             --  x-coordinate of the vertex
                - y (float)             --  y-coordiante of the vertex
                - connections (list)    --  list of coordiante tuples that are
                                            connected to this vertex
        '''
        self.x = vert_obj.x
        self.y = vert_obj.y
        self.connections = vert_obj.connections
        self.rem_req = 0

        # If a vertex has no connections, it's not useful for defining a line
        if len(self.connections) == 0:
            msg = '''The vertex, {}, to be activated does not have any connections
                     and should be removed.'''.format(vert_obj)
            raise RuntimeError(msg)

    def __repr__(self):
        return 'Active Vertex({}, {})'.format(self.x, self.y)

    def request_removal(self):
        '''
        Request that this vertex be deactivated. When the number of remove
        requests equals the number of connections to a vertex, the vertex will
        raise a ValueError.

        RAISES:
        ValueErorr  --  When the number of remove requests equals or exceeeds
                        the number of connections for the vertex.
        '''
        self.rem_req += 1 # Add a remove request
        if self.rem_req >= len(self.connections):
            msg = '''Remove requests, {}, exceed number of connections, {}. This
                     vertex should be deactivated
                     '''.format(self.rem_req, len(self.connections))
            raise ValueError(msg)

class active_vert_dict(SortedDict):
    '''
    Dictionary of active vertices that also contains a method for effectively
    issuing remove requests and removing the vertex from the dictionary and
    sort primarily by the y-coordinate (reverse of default behavior)

    Special attribute:
    split_lines (dict)      --  Active connected vertices that have been split
                                in the hierarchy by the introduction of a new
                                vertex. These are special candidates for further
                                intersection testing since interesctions could
                                go unnoticed if the end points are far away from
                                the intersection area. The dict is structured
                                such that KEY = end-point coordinates and
                                VALUE = other end-point coordinates. In this
                                way, the line is listed twice, once from each
                                direction.
    '''
    def __init__(self, *args, **kwargs):
        '''
        Initializes as a SortedDict designed to sort by y-coordinate and then
        by x-coordinate (reverse of default behavior)
        '''
        SortedDict.__init__(self, lambda tup: (tup[1],tup[0]), *args, **kwargs)
        self.split_lines = {}

    def __repr__(self):
        return SortedDict.__repr__()

    def request_removal(self, vertex_coords):
        '''
        Create a request that the vertex at vertex_coords be deactivated. The
        vertex will only be deactiviated once the removal requests equal the
        number of connections for the vertex.
        '''
        try:
            self[vertex_coords].request_removal()
        except ValueError:
            # Number of removal requests equal connections
            if vertex_coords in self.split_lines:
                other_end_coords = split_lines[vertex_coords]
                self.split_lines.pop(vertex_coords)
                self.split_lines.pop(other_end_coords)
                # Both points are removed in this step even though one would
                #  assume the other point would be subsequently removed. However,
                #  it seems more important to not break the split lines dict by
                #  having a dangling point somehow. Perhaps this is unnecessary.
            return {vertex_coords:self.pop(vertex_coords)}
        except KeyError:
            raise KeyError('{} is not an active vertex'.format(vertex_coords))

    def split_verts(self, v_coords1, v_coords2):
        '''
        Identify that the line formed between v_coords1 and v_coords2 is split
        by at least one other vertex.
        '''
        self.split_lines.update({v_coords1:v_coords2, v_coords2:v_coords1})

def distance(point1, point2):
    '''Calculates the distance between two points'''
    dist = math.sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2)
    return dist

def check_ccw(A, B, C, tol=1e-06):
    '''
    Checks whether three points A, B, and C are oriented counter-clockwise
    relative to each other. There are three possible outcomes:
     1  - The points form a counter-clockwise (CCW) triangle
    -1  - The points form a clockwise (CW) triangle
     0  - The points are colinear within the specified tolerance, tol

     Note: the tolerance is the maximum x and y distance a point could be
     away from a line while still being considered colinear. Implicitly, this
     means the maximum distance a point could be from a line would be
     sqrt(2)*tol when the line AB is at a +/-45 degree angle to the horizontal.

     ARGUMENTS:
     A, B, C (iter)     --  three points each given by a length-2 iterable

     EXCEPTIONS:
     ValueError         --  if all three points are the same
    '''
    # Check that points differ
    if (A[0] == B[0]) and (A[1] == B[1]) and (A[0] == C[0]) and (A[1] == C[1]):
        raise ValueError('{}, {}, and {} have the same coordiantes'.format(A,B,C))

    # Create test value (a difference of slopes)
    test = (C[1] - A[1])*(B[0] - A[0]) - (B[1] - A[1])*(C[0] - A[0])

    # Check for colinearity within worst-case position (NOT slope) tolerance 
    if abs(test) <= tol*(B[0] - A[0] + B[1] - A[1]):
        return 0
    # Check for counter-clockwise
    elif test >= 0:
        return 1
    # Otherwise it's clockwise
    else:
        return -1

def intersection(line1, line2, tol=1.0e-06):
    '''
    Checks whether two lines intersect each other by comparing the chirality of
    all combinations of three points taken from the four end points. Five
    possible cases are idenitifed:
    1) the lines form a simple intersection with no three points being colinear
    2) two lines share a common end-point
    3) a point on one line is colinear with the points on the other line
    4) both lines are colinear
    5) the lines do not intersect at all

    If the lines interesect, the number of the case (1-4) will be returned
    along with the results of the counter-clockwise tests that help idenitfy
    how the two lines intersect. If the lines do not intersect, nothing is
    returned.

    ARGUMENTS:
    line1 (tuple 2) --  length-2 tuple containing the end points of the first
                        line with each point each defined as an ordered pair of
                        coordinates
    line2 (tuple 2) --  the second line with the same syntax as the first line

    OPTIONAL ARGUMENTS:
    tol (float)     --  tolerance to which any three points will be evaluated as
                        either counter-clockwise or clockwise. In this sense, a
                        point can be a maximum distance of sqrt(2)*tol from a
                        line before being considered not coincident on that line

    RETURNS:
    None            --  If the lines do not intersect
    case (int)      --  The case number above that is identified for the two
                        lines except when the lines do not intersect
    ccw_tests (list)--  List of the individual tests for chirality. For example,
                        if two lines are defined as AB and CD, the list will be
                        results for [ABC, ABD, CDA, CDB] chirality. If the
                        the points are arragned like {ABC: ccw, ABD: cw, 
                        CDA: cc, CDB: ccw}, the list will be [1, -1, -1, 1] and
                        it is apparent that the lines intersect. If any value is
                        zero, that means that those three points were found to
                        be colinear.

    EXCEPTIONS:
    ValueError      --  if lines are the same or reversed
    '''
    # Lines are referred to as AB and BC
    a = line1[0]
    b = line1[1]
    c = line2[0]
    d = line2[1]

    # Check to make sure lines are not the same or reversed
    if (a == c and b == d) or (b == c and a == d):
       raise ValueError('Lines {} and {} are colinear'.format(line1, line2))

    # Test three points at a time for chirality
    ccw_tests = [check_ccw(a, b, c, tol=tol),
                 check_ccw(a, b, d, tol=tol),
                 check_ccw(c, d, a, tol=tol),
                 check_ccw(c, d, b, tol=tol)]

    # Sum absolute values of tests to determine how many points are colinear
    ccw_sum = sum([abs(x) for x in ccw_tests])

    # Check for the type of intersection
    if ccw_sum == 0:
        # All Points are colinear
        points = {'a':a, 'b':b, 'c':c, 'd':d}
        sort_points = sorted(points, key=points.__getitem__)
        pairs = {'a':'b', 'b':'a', 'c':'d', 'd':'c'}
        #Check if paired points are adjacent after being sorted
        if sort_points[1] != pairs[sort_points[0]]:
            case = 4 # All points are colinear and lines overlap
        else:
            return None # Points are colinear but do not overlap
    elif ccw_sum == 2:
        case = 2 # An endpoint is shared between lines
    elif ccw_sum == 3:
        # Three points are colinear but don't necessarily intersect
        if sum(ccw_tests[0:2]) == 0 or sum(ccw_tests[2:4]) == 0:
            # An intersection occurs when abc and abd alternate chirality OR
            #  when cda and cdb alternate chirality
            case = 3
        else:
            # Otherwise three points are colinear but no intersection
            return None
    else:
        # No three points are colinear
        if sum(ccw_tests[0:2]) == 0 and sum(ccw_tests[2:4]) == 0:
            # An intersection occurs when abc and abd alternate chirality AND
            #  when cda and cdb alternate chirality
            case = 1
        else:
            # Otherwise three points are colinear but no intersection
            return None
    # Return the situation
    return case, ccw_tests

def find_intersections(vertices, tol=1e-06):
    '''
    Finds the intersections between line segments. These intersections can
    be of four forms:
    1) two lines simply intersect and have different slopes
    2) two lines share a common end-point
    3) one endpoint lies on the other line
    4) both lines are colinear and overlap

    or no intersection occurs.

    A sweep-line algorithm is used to locate
    intersections. This works by first ordering the vertices from left to right/
    bottom to top and then looping through the vertices. A vertex is "activated"
    when the sweep line passes through it. If the vertex is a left-bottom
    vertex, it is added into the list of active vertices for the first time.
    Then the connections for this vertex are also added to the list of active
    vertices. The algorithm then checks for intersections with all lines
    connected to vertices that are between the line endpoints, vertices that are
    above or below the current line, or any lines that are still active and
    have previously been split by other vertices.
    If a vertex is already active when it is reached by the sweep-line, any
    inactive connections are added to the active hierarchy and remove requests
    are issued for the vertex and its active connections. When the number of
    remove requests equals the number of connections on a vertex, the vertex
    is deactivated, meaning that the sweep-line has already passed over all
    vertices connected to the deactivated vertex.

    The intersections are stored as a tuple of the two lines that have
    intersected in addition to the type of intersection that has ocurred.

    ARGUMENTS:
    vertices (dict)         --  Dictionary (or object that can be converted into
                                a dictionary) of vertex coordinates and vertex
                                objects that satisfy the requirements necessary
                                to be converted into active vertex objects.

    OPTIONAL ARGUMENTS:
    tol (float)             --  Tolerance to which intersection algorithm is run
                                (Default: 1e-06)

    RAISES:
    RuntimeError            --  for the following fatal errors:
                                - The vertex coordinates are inconsistent with
                                  the x and y attributes of the vertex object
                                - Connection information between vertices is
                                  inconsistent. In other words, one vertex says
                                  it's connected to another vertex but that
                                  other vertex doesn't agree.
                                - A new vertex has connections that are already
                                  active. This indicates that the algorithm is
                                  flawed because if a vertex is new, then by
                                  definition it shouldn't be connected to any
                                  vertices that were previously active.
    '''
    # First order the verticies from left to right (i.e by x coordinate)
    sorted_verts = SortedDict(vertices)

    # Initialization
    active_verts = active_vert_dict()
    removed_vertices = {}
    prev_vertex_pos = None
    intersection_list = []

    # Sweep line algorithm from left to right
    for vertex_coords, vertex in sorted_verts.iteritems():

        # Check for consistency
        if vertex_coords != (vertex.x, vertex.y):
            msg = '''FATAL ERROR: vertex coordinates refer to wrong vertex\n
                     coordinates: {}\t (x,y): ({},{})'''.format(vertex_coords, 
                     vertex.x, vertex.y)
            raise RuntimeError(msg)

        if prev_vertex_pos: # Check whether it's the first iteration
            if vertex.x - tol >= prev_vertex_pos[0]:
                removed_vertices.clear()

        # Update the the position of the sweep line
        prev_vertex_pos = vertex_coords

        # If the vertex is already active, it is a right vertex
        if vertex_coords in active_verts:
            # Use a simpler notation
            existing_vert = active_verts[vertex_coords]

            # Loop through connections
            for con in existing_vert.connections:
                if con in active_verts:
                    # Request removal for active connection
                    removed_vertices.update(active_verts.request_removal(con))
                    # Make removal request on current vertex for every active
                    #  connection
                    removed_vertices.update(active_verts.request_removal(vertex_coords))
                else:
                    # Otherwise add connection to active vertex list
                    active_verts.update({con:active_vertex(sorted_verts[con])})

        # Since the vertex is new, it is a left vertex
        else:
            # Begin by adding the vertex to the active dict
            active_verts.update({vertex_coords:active_vertex(vertex)})

            # Did adding this vertex split an existing line?
            loc = active_verts.index(vertex_coords)
            pred_coord = active_verts.iloc[loc-1]
            succ_coord = active_verts.iloc[loc+1]

            # Check whether new vertex splits two connected active vertices
            if pred_coord in active_verts[succ_coord].connections:
                # Consistency check
                if succ_coord not in active_verts[pred_coord].connections:
                    msg = '''FATAL ERROR: vertex {} shows connection with {}
                             but not visa versa'''.format(succ_coord, pred_coord)
                    raise RuntimeError(msg)
                # Line has been split, add to dict of split lines
                active_verts.split_verts(pred_coord, succ_cord)

            # Now loop through the connections of the new vertex (all of which
            #  themselves must be new)
            for con in vertex.connections:
                # Consistency check
                if con in active_verts:
                    msg = '''FATAL ERROR: Vertex connections for new vertices
                             should all be inactive. If there are any active
                             vertex connections, then the vertex shouldn\'t be
                             new. New vertex: {}, active connection: {}
                             '''.format(vertex_coords, con)
                    raise RuntimeError(msg)

                # Add the connection as an active vertex
                active_verts.update({con:active_vertex(sorted_verts[con])})

                # Index the vertex and connection
                vert_loc = active.verts.index(vertex_coords)
                con_loc = active_verts.index(con)

                # Identify predecessor and successor indices
                if vert_loc < con_loc:
                    pred_loc = vert_loc - 1
                    succ_loc = con_loc + 1
                else:
                    pred_loc = con_loc - 1
                    succ_loc = vert_loc + 1

                # Is this connection already split?
                if abs(con_loc - vert_loc) > 1:
                    # Yes this line is split upon being added
                    active_verts.split_verts(vert_loc, con_loc)
                    # Loop through vertices lying between line endpoints
                    internal_points = active_verts.islice(pred_loc + 2, succ_loc - 1)
                    for i_pt_coords in internal_points:
                        # Loop through connections on each internal vertex
                        for i_pt_con in active_verts[i_pt_coords].connections:
                            # Filter by active vertices
                            if i_pt_con in active_verts:
                                # Check intersection between current line and
                                # this connection on the internal point
                                line1 = (vertex_coords, con)
                                line2 = (i_pt_coords, i_pt_con)
                                result = intersection(line1, line2, tol)
                                # Will not append if None is returned
                                if result:
                                    intersection_list.append((line1, line2), result)

                # Now look at the predecessor to the newly inserted line segment
                pred_coords = active_verts.iloc[pred_loc]
                # Loop through connections
                for pred_con in active_verts[pred_coords].connections:
                    line1 = (vertex_coords, con)
                    line2 = (pred_coords, pred_con)
                    result = intersection(line1, line2, tol)
                    if result:
                        intersection_list.append((line1, line2), result)

                # Find intersections with the successor now
                succ_coords = active_verts.iloc[suc_loc]
                for succ_con in active_verts[succ_coords].connections:
                    line1 = (vertex_coords, con)
                    line2 = (succ_coords, succ_con)
                    result = intersection(line1, line2, tol)
                    if result:
                        intersection_list.append((line1, line2), result)

                # Look through connections of removed vertices and find
                #  intersections
                for rem_vert_coords, rem_vertex in removed_vertices.iteritems():
                    for rem_vert_con in rem_vertex.connections:
                        line1 = (vertex_coords, con)
                        line2 = (rem_vert_coords, rem_vert_con)
                        result = intersection(line1, line2, tol)
                        if result:
                            intersection_list.append((line1, line2), result)

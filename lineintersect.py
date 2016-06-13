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
                     and should be removed.'''.format{vert_obj}
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
            # when cda and cdb alternate chirality
            case = 3
        else:
            # Otherwise three points are colinear but no intersection
            return None
    else:
        # No three points are colinear
        if sum(ccw_tests[0:2]) == 0 and sum(ccw_tests[2:4]) == 0:
            # An intersection occurs when abc and abd alternate chirality AND
            # when cda and cdb alternate chirality
            case = 1
        else:
            # Otherwise three points are colinear but no intersection
            return None
    # Return the situation
    return case, ccw_tests

def find_intersections(vertices):
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
    when the sweep line passes through it. If the vertex is a left or bottom
    vertex, it is added into the list of active vertices. If a vertex is a right
    or top vertex, the algorithm checks for an intersection between all active
    vertices that vertically lie between the two endpoints plus one vertex above
    the line and one below. If the right/top point has other connections, it is
    added to the active vertex list and a remove request is made against the
    left/bottom vertex. Once the number of remove requests equals the number of
    connections on that vertex, that vertex is deactivated.

    The intersections are stored as a tuple of the two lines that have
    intersected in addition to the type of intersection that has ocurred.
    '''
    # First order the verticies from left to right (i.e by x coordinate)
    sorted_verts = SortedDict(vertices)






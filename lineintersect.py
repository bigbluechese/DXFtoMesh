#!/usr/bin/env python

'''
Calculates the intersections points of a bunch of lines and then breaks those
lines up into segments such that no lines intersect and instead only share
end points
'''

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
    all combinations of three points taken from the four end points. Three
    possible cases are idenitifed:
    1) the lines form a simple intersection with no three points being colinear
    2) two lines share a common end-point
    3) a point on one line is colinear with the points on the other line
    4) both lines are colinear
    5) the lines do not intersect at all

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

    # Check for the type
    if ccw_sum == 0:
        case = 4 # All points are colinear
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

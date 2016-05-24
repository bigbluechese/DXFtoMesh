#!/usr/bin/env python

'''
Collection of helper functions used with DXftoSegments module
'''

from __future__ import print_function
import numpy as np
import re
from numbers import Number

def angle360(angle, decimals=15):
    '''Converts an angle from a domain of [-pi, pi] to [0, 2*pi]'''
    angle = np.around(angle, decimals)
    if angle < 0:
        new_angle = 2*np.pi + angle
    elif angle >= 2*np.pi:
        new_angle = angle - 2*np.pi
    else:
        new_angle = angle
    # Correct for negative zero angles
    if new_angle == -0.0:
        new_angle = 0.0
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

def approx(val, tol=1.0e-08):
    '''Approximates a value to the given tolerance'''
    decimals = np.abs(np.floor(np.log10(np.abs(tol))).astype(int))
    new_val = round(val, decimals)
    if new_val == -0.0:
        new_val = 0.0 #Get rid of negative zero
    #print val, new_val
    return new_val


def bulge_to_arc(start, end, bulge, tol=1.0e-08):
    '''
    Converts bulge information into arc information and outputs all of the
    information as a tuple.

    ARGUMENTS:
    start (tuple)       --  Starting coordinates for segment
    end (tuple)         --  Ending coordinates for segment
    bulge (float)       --  Bulge value associated with the segment as defined
                            by Autodesk conventions

    RETURNS:
    ((start, end), (bulge, start_angle, end_angle, center, radius))
    start (tuple)       --  See arguments
    end (tuple)         --  See arguments
    bulge (float)       --  See arguments (negative for clockwise arc)
    start_angle (float) --  Starting angle for arc
    end_angle (float)   --  Ending angle for arc
    center (tuple)      --  Coordinates for center of arc
    radius (float)      --  Radius of arc
    '''
    # Distance between points
    d = np.sqrt((start[0] - end[0])**2 + (start[1] - end[1])**2)
    # Angle between points from center
    theta = 4*np.arctan(bulge)
    # Radius of circle making arc
    radius = approx(d/2/np.sin(abs(theta)/2), tol=1e-12)
    # Find angle of segment relative to x axis
    alpha = np.arctan2(end[1]-start[1], end[0]-start[0])
    # Find angle between radius vector (from center to point 1) and x-axis
    if bulge > 0:
        # Find angle between segment and radius. Beta is negative if theta is
        # greater than pi
        beta = np.pi/2 - theta/2
        # Angle to radius vector from x-axis is then the SUM of alpha and beta
        gamma = alpha + beta
    else:
        # Find angle between segment and radius. Beta is negative if theta is
        # greater than pi
        beta = np.pi/2 + theta/2
        # Angle to radius vector from x-axis is then the DIFFERENCE between
        # alpha and beta
        gamma = alpha - beta
    # Gamma angle and radius describe the vector pointing from the start point
    # to the center
    center = (approx(radius*np.cos(gamma)+start[0], tol=tol),
              approx(radius*np.sin(gamma)+start[1], tol=tol))
    # Now compute start and stop angles relative to horizontal in a
    # counter-clockwise sense
    start_angle = angle360(np.arctan2(start[1]-center[1], 
                                      start[0]-center[0]))
    end_angle = angle360(np.arctan2(end[1]-center[1],
                                    end[0]-center[0]))

    # Compile all bulge/arc information and add it to segments
    return ((start, end), (bulge, start_angle, end_angle, center, radius))

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

# def tuple_string_check(var):
#     '''
#     Checks whether a given input is a str conversion of a 2-length tuple and
#     then ensures that the numbers in the tuple are printed like floats. If a 
#     tuple is given instead of a string, it will convert the tuple to the 
#     correct string format.

#     ARGUMENTS:
#     var (str)               --  variable to be tested

#     RETURNS:
#     str_tup (str)           --  properly formatted string of tuple

#     RAISES:
#     TypeError               --  if var is not a string or tuple
#     ValueError              --  if var is not a length-2 tuple that contained
#                                 numbers
#     '''
#     if type(var) == tuple:
#         var = str((float(var[0]), float(var[1])))
#     elif type(var) != str:
#         raise TypeError('variable must be a string conversion of a tuple of two numbers')
    
#     # Compile regular expression and check whether tuple is length 2 containing
#     # numbers
#     tuple_check = re.compile('[(]([-+]?[0-9]*\.?[0-9]+[eE]?[-+]?[0-9]*), ([-+]?[0-9]*\.?[0-9]+[eE]?[-+]?[0-9]*)[)]')
#     tuple_match = tuple_check.match(var)
#     if tuple_match:
#         tuple_float = (float(tuple_match.groups()[0]), float(tuple_match.groups()[1]))
#     else:
#         raise ValueError('tuple converted to string must have been length 2 and had only numbers')
    
#     # Return the new string that is a tuple of two floats
#     str_tup = str(tuple_float)
#     return str_tup

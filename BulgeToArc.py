#!/usr/bin/python

import math
import numpy as np

def distance(segment):
    start = segment[0][0]
    end = segment[0][1]
    d = math.sqrt((start[0] - end[0])**2 + (start[1] - end[1])**2)
    return d

def bulge_to_angle(segment):
    bulge = segment[1][0]
    angle = 4*math.atan(bulge)
    if angle < 0:
        new_angle = angle + 2*math.pi
    else:
        new_angle = angle
    print 'raw angle: {} new angle: {}°'.format(math.degrees(angle), math.degrees(new_angle))
    return new_angle

def radius(segment):
    bulge = segment[1][0]
    angle = 4*math.atan(bulge)
    d = distance(segment)
    radius = d/2/math.sin(abs(angle)/2)
    return radius

def beta(segment):
    bulge = segment[1][0]
    angle = 4*math.atan(bulge)
    beta = (math.pi/2 - abs(angle)/2)*angle/abs(angle)
    return math.degrees(beta)

def beta2(segment):
    bulge = segment[1][0]
    theta = 4*math.atan(bulge)
    beta = -(math.pi/2 - abs(theta)/2)*(math.pi - theta)/abs(math.pi - theta)   
    return math.degrees(beta)

def center(segment):
    start = segment[0][0]
    end = segment[0][1]
    r = radius(segment)
    alpha = math.atan2(end[1]-start[1], end[0]-start[0])
    beta = math.radians(beta2(segment))
    gamma = alpha - beta
    print math.degrees(gamma)
    center = (r*math.cos(gamma)+start[0], r*math.sin(gamma)+start[1])
    return center



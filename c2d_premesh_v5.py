from numpy import array, mean, linalg, arccos, vdot, isnan, pi, cross, random, ones
import time
from math import atan2
import numpy as np
from matplotlib.collections import PolyCollection
import matplotlib as mpl
#mpl.use('TkAgg')
import matplotlib.pyplot as plt
import pdb
import dxfgrabber
import pickle
import DXFtoSegments

'''
Jeff's notes:
Interactive plotting isn't working on Mac OSX. Also the mpl.use function was
placed after the pyplot import which apparently did nothing.

Default plotting behavior has been changed from True to False
'''


class Vertex(object):

    def __init__(self, id, coords):
        self.id = id
        self.coords = coords
        self.region_assoc = set()
        self.max_connections = []
        self.connections = set()
        self.boundary = False

        
    def __str__(self):
        return "Vertex " + str(self.id) + " has coordinates " + str(self.coords)
        
    def __repr__(self):
        return "{}_{}".format(self.__class__.__name__, self.id)
        
    def associateRegion(self, regionID, corner):
        self.region_assoc.append([regionID, corner])
        
    

class Edge(object):
    
    def __init__(self, id, Vertex0, Vertex1):
        self.id = id
        self.Vertex = [Vertex0, Vertex1]
        self.region_assoc = []
        self.boundary = False

    
    def __str__(self):
        return "Edge " + str(self.id) + " has Vertices " + str(self.Vertex[0].id) + " and " + \
        str(self.Vertex[1].id)
        
    def __repr__(self):
        return "{}_{}".format(self.__class__.__name__, self.id)
        
    def associateRegion(self, regionID, side, direction):
        self.region_assoc.append([regionID, side, direction])
        
    def unit_vec(self, starting_vertex):
        if starting_vertex == self.Vertex[0]:
            return (self.Vertex[1].coords - self.Vertex[0].coords)/\
            linalg.norm(self.Vertex[1].coords - self.Vertex[0].coords)
        elif starting_vertex == self.Vertex[1]:
            return (self.Vertex[0].coords - self.Vertex[1].coords)/\
            linalg.norm(self.Vertex[0].coords - self.Vertex[1].coords)
            
    def other_vertex(self, one_vertex):
        if one_vertex == self.Vertex[0]:
            return self.Vertex[1]
        elif one_vertex == self.Vertex[1]:
            return self.Vertex[0]
            
    def is_a_vertex(self, one_vertex):
        if one_vertex == self.Vertex[0]:
            return 0
        elif one_vertex == self.Vertex[1]:
            return 1
        else:
            return -1
            
    @property
    def mdpt(self):
        return ((self.Vertex[0].coords[0]+ self.Vertex[1].coords[0])/2.0, \
        (self.Vertex[0].coords[1] + self.Vertex[1].coords[1])/2.0)
    
    @property    
    def max_regions(self):
        if self.boundary:
            return 1
        else:
            return 2

class Region(object):
    
    unit_vecs = [array([0,1]), array([1,0]), array([0,-1]), array([-1,0])]
    
    def __init__(self, id, list_of_Edges):
        self.id = id
        self.Edges={}
        #self.Vertices={}

        # The meshregion mdpt is the average position of all its vertices
        self.mdpt = (mean([edge.mdpt[0] for edge in list_of_Edges]), \
        mean([edge.mdpt[1] for edge in list_of_Edges]))
        
        # Find vectors from meshregion mdpt to each of the edge mdpts
        self.edge_mdpts = [edge.mdpt for edge in list_of_Edges]
        self.edge_mdpt_vecs = [(array(edge.mdpt) - array(self.mdpt))/linalg.norm((array(edge.mdpt) - array(self.mdpt))) for edge in list_of_Edges]
    
        # For each cardinal direction (N,E,S,W)
        for i in range(len(self.unit_vecs)):
        
            # Find which (meshregion mdpt)-to-(edge-mdpt) vector has the smallest error from
            # the cardinal direction in question by taking the norm of the difference of the 
            # two vectors
            self.cardinal_norms = array([linalg.norm(self.unit_vecs[i] - edge_mdpt_vec) for edge_mdpt_vec in self.edge_mdpt_vecs])
    
            if i == 0: # North
                self.Edges['n'] = list_of_Edges[self.cardinal_norms.argmin()]
                # If Vertex[0]'s x-coord is greater than Vertex[1]'s x-coord
                if self.Edges['n'].Vertex[0].coords[0] > self.Edges['n'].Vertex[1].coords[0]:
                    self.Edges['n'].associateRegion(self.id, 'n', 1)
                else:
                    self.Edges['n'].associateRegion(self.id, 'n', -1)
                    
            elif i == 1: # East
                self.Edges['e'] = list_of_Edges[self.cardinal_norms.argmin()]
                # If Vertex[1]'s y-coord is greater than Vertex[0]'s y-coord
                if self.Edges['e'].Vertex[1].coords[1] > self.Edges['e'].Vertex[0].coords[1]:
                    self.Edges['e'].associateRegion(self.id, 'e', 1)
                else:
                    self.Edges['e'].associateRegion(self.id, 'e', -1) 
                          
            elif i == 2: # South
                self.Edges['s'] = list_of_Edges[self.cardinal_norms.argmin()]
                # If Vertex[1]'s x-coord is greater than Vertex[0]'s x-coord
                if self.Edges['s'].Vertex[1].coords[0] > self.Edges['s'].Vertex[0].coords[0]:
                    self.Edges['s'].associateRegion(self.id, 's', 1)
                else:
                    self.Edges['s'].associateRegion(self.id, 's', -1) 
                    
            elif i == 3: # West#
                self.Edges['w'] = list_of_Edges[self.cardinal_norms.argmin()]
                # If Vertex[0]'s y-coord is greater than Vertex[1]'s y-coord
                if self.Edges['w'].Vertex[0].coords[1] > self.Edges['w'].Vertex[1].coords[1]:
                    self.Edges['w'].associateRegion(self.id, 'w', 1)
                else:
                    self.Edges['w'].associateRegion(self.id, 'w', -1) 
                    
                
        # Find matching vertices in adjacent edges to determine vertex directions (ne, nw, se, sw)  
        
        #self.Vertices['ne'] = set(self.Edges['n'].Vertex).intersection(self.Edges['e'].Vertex).pop()
        #self.Vertices['nw'] = set(self.Edges['n'].Vertex).intersection(self.Edges['w'].Vertex).pop()
        #self.Vertices['se'] = set(self.Edges['s'].Vertex).intersection(self.Edges['e'].Vertex).pop()
        #self.Vertices['sw'] = set(self.Edges['s'].Vertex).intersection(self.Edges['w'].Vertex).pop()
        
    def Vertices(self, ord_dir):
        try:
            return set(self.Edges[ord_dir[0]].Vertex).intersection(self.Edges[ord_dir[1]].Vertex).pop()
        except Exception:
            print "Not a valid direction"
            return 0
                        
        
    def Neighbors(self, card_dir):
        # Find neighbors
        try:
            return [region_info for region_info in self.Edges[card_dir].region_assoc if region_info[0] != self.id][0]
        except IndexError:
            #print "Region " + str(self.id) + " has no neighbor in the " + card_dir + " direction."
            return 0
        except KeyError:
            #print card_dir + " is not a valid direction"
            return None

def angle_between(v1, v2):
    """ Returns the clockwise angle in radians (0 to 2pi) between vectors 'v1' and 'v2':: """
    
    consang = vdot(v1, v2)
    sinang = cross(v1, v2)
    angle = atan2(sinang, consang)
    
    if angle < 0:
        angle += 2*pi
    return angle
         
def find_most_CCW(starting_vertex, leading_in_Edge, Edge_list):
    angle_list = array([angle_between(edge.unit_vec(starting_vertex), leading_in_Edge.unit_vec(starting_vertex)) for edge in Edge_list])
    
    return Edge_list[angle_list.argmin()]
    
def find_most_CW(starting_vertex, leading_in_Edge, Edge_list):
    angle_list = array([angle_between(edge.unit_vec(starting_vertex), leading_in_Edge.unit_vec(starting_vertex)) for edge in Edge_list])
    
    return Edge_list[angle_list.argmax()]

def find_most_SW(Vertex_list):
    vertex_array = array([list(vertex.coords) for vertex in Vertex_list])
    
    # The SW bounding box corner is the point with minimum x and minimum y value of any vertex
    SW_bounding_box_corner = [min(array(vertex_array)[:,0]), min(array(vertex_array)[:,1])]
    
    a = SW_bounding_box_corner[0]*ones(len(vertex_array))
    b = SW_bounding_box_corner[1]*ones(len(vertex_array))
    
    # Compute the vector from SW BB corner from each vertex
    vector_from_SW_BB = vertex_array - array([a, b]).T
    
    # The vertex whose distance is closest to the SW BB corner is considered the Southwestern-most point
    # and is returned
    norms = array([linalg.norm(vertex - SW_bounding_box_corner) for vertex in vertex_array])
    return Vertex_list[norms.argmin()]
    
class C2DMesh(object):

    
    def __init__(self, vertex_list, edge_list, plotting=False, plotting_pause=0.5):
        # Preallocate space for a list of Vertex objects   
        self.Vertex_list = [None]*len(vertex_list)
        self.Edge_list = [None]*len(edge_list)
        self.plotting_pause = plotting_pause
        
        self.plotting = plotting
        
        # We want some code here to check that the edges make sense, like they only use vertices specified
        #
        #
        
        for i in range(len(vertex_list)):
            self.Vertex_list[i] = Vertex(i, array(vertex_list[i]))
        
        for i in range(len(edge_list)):
            self.Edge_list[i] = Edge(i, self.Vertex_list[edge_list[i][0]], self.Vertex_list[edge_list[i][1]])
            
        # Fill in each vertex's .max_connections
        for i in range(len(vertex_list)):
            self.Vertex_list[i].max_connections = len([edge for edge in self.Edge_list if edge.is_a_vertex(self.Vertex_list[i]) >=0 ])
            print 'Vertex ' + str(i) + ' has ' + str(self.Vertex_list[i].max_connections) + ' connections.'
                  
        current_mesh_id = 1
        edges_in_current_region = []
        self.Region_list = []
        active_edge = None
        active_vertex = None
        
        outer_loop_i = 0
        inner_loop_i = 0
        restart_search = 0
        
        # Plotting
        if self.plotting:
            plt.close('all')
            fig, ax = plt.subplots()
            plt.ion()
            plt.show()
        
        
        self.mark_boundary_edges()
           
        #pdb.set_trace()     
        # Outer loop to find mesh regions.  The outer loop determines which vertex and edge to first look at
        # when trying to find a new mesh region.  For the first mesh region, it just finds __.  For subsequent
        # mesh regions, it builds off edges that are already associated with a previously found mesh region.        
        
        while True:
            print "Starting search for mesh region {}".format(current_mesh_id)
            if len(self.Region_list) == 0:
                # Let's try to start at the Southwestern-most vertex.
                active_vertex = find_most_SW(self.Vertex_list)
                print "First vertex: ", active_vertex
                
                next_edge_candidates = [edge for edge in self.Edge_list if edge.is_a_vertex(active_vertex) >=0 and edge != active_edge]
                
                temp_edge = Edge(-1, active_vertex, Vertex(-1, active_vertex.coords + array([.01,0])))
                active_edge = find_most_CCW(active_vertex, temp_edge, next_edge_candidates)
                
                #edges_in_current_region.append(active_edge)
                
            else:
                # Query vertices that have not maxed out their max_connections
                available_vertices = [vertex for vertex in self.Vertex_list if len(vertex.connections) < vertex.max_connections]

                if len(available_vertices) == 0:
                    print "All possible vertices have been used!"
                    break
                
                # Edge search 1: find edges that are only associated with one region, so far.
                edge_search1 = [edge for edge in self.Edge_list if len(edge.region_assoc) == 1]
                edge_search2 = [edge for edge in edge_search1 if len(edge.region_assoc) < edge.max_regions]
                active_edge = edge_search2[restart_search]
                
                # If the first region this edge is associated with is running in one direction, then the adjacent region who shares that edge
                # should run in the opposite direction.  See: Winged-edge
                if active_edge.region_assoc[0][2] > 0:
                    active_vertex = active_edge.Vertex[1]
                else:
                    active_vertex = active_edge.Vertex[0]
                
            print "First vertex: ", active_vertex
        
            inner_loop_i = 0
        
            while True:
                
                # The new active vertex is the other end of the current edge
                active_vertex = active_edge.other_vertex(active_vertex)
                print active_vertex
                
                # The next edge candidates are ones that are associated with the new active vertex, with the
                # exception of the current active edge.
                next_edge_candidates = [edge for edge in self.Edge_list if edge.is_a_vertex(active_vertex)>= 0\
                and edge != active_edge]
                #print next_edge_candidates
                
                # The next active edge is the one with the smallest clockwise angle, or most counter-clockwise
                active_edge = find_most_CCW(active_vertex, active_edge, next_edge_candidates)
                
                # Append the new active edge to the list of edges in the current mesh region
                edges_in_current_region.append(active_edge)
                #print edges_in_current_region
                
                inner_loop_i +=1
                #print inner_loop_i
            
                # If we've found four edges, we found a region!
                if len(edges_in_current_region) >= 4:
                    
                    # Need some code here to make sure we have indeed found a closed quadrilateral
                    if len(set( [edge.Vertex[i] for i in range(2) for edge in edges_in_current_region])) != 4:
                        print "THIS IS NOT A QUADRILATERAL"
                        time.sleep(1)
                        restart_search += 1
                        # Clear edges_in_current region
                        edges_in_current_region = []
                        #raise RuntimeError('This is not a quadrilateral')
                        
                    # Create a meshregion object and add it to Region_list
                    m = Region(current_mesh_id, edges_in_current_region)
                    self.Region_list.append(m)
                    
                    # Associate each vertex with edges that are connected to it
                    # There is a better way to do this than to iterate through each vertex
                    for i in range(len(self.Vertex_list)):
                        for k in range(len(edges_in_current_region)):
                            if edges_in_current_region[k].is_a_vertex(self.Vertex_list[i])>=0:
                                self.Vertex_list[i].connections.add(edges_in_current_region[k].other_vertex(self.Vertex_list[i]).id)
                                self.Vertex_list[i].region_assoc.add(m.id)
                            
                    # Clear edges_in_current region
                    edges_in_current_region = []
                    
                    # Plotting
                    if self.plotting:
                        plot_vertices = []
                    
                        plot_vertices.append([m.Vertices('ne').coords,
                                              m.Vertices('nw').coords,
                                              m.Vertices('sw').coords,
                                              m.Vertices('se').coords])
                        
                        plt.annotate(str(m.id), xy=(m.mdpt[0], m.mdpt[1]))
                        
                        plot_vertices = array(plot_vertices)
                        
                        #z = array(range(len(self.Region_list)))
                        z = array(current_mesh_id)
                                
                        #Make the collection and add it to the plot.
                        coll = PolyCollection(plot_vertices, array=z, edgecolors='000000')
                        coll = PolyCollection(plot_vertices, array=z, cmap=mpl.cm.Pastel1, edgecolors='000000')
                        ax.add_collection(coll)
                        ax.autoscale_view()
                        ax.set_aspect('equal', 'datalim')
                        plt.draw()
                        plt.pause(self.plotting_pause)
                        plt.show()
            
                    print 'Mesh region {} found!'.format(current_mesh_id)

                    current_mesh_id += 1 # increment the current region id number
                    restart_search = 0
                    break
                    
                
        print "\n\n\n Done finding regions! \n"
        #print self.Region_list[0].Edges.values()
        #print self.Region_list[1].Edges.values()
        #print self.Region_list[2].Edges.values()

    def mark_boundary_edges(self):
        active_edge = None
        active_vertex = None
        
        # Let's try to start at the Southwestern-most vertex.
        boundary_finder_i = 0
        while True:
            if boundary_finder_i == 0:
                active_vertex = find_most_SW(self.Vertex_list)
                active_vertex.boundary = True       
                #print "active vertex: ", active_vertex
            
                next_edge_candidates = [edge for edge in self.Edge_list if edge.is_a_vertex(active_vertex) >=0 and edge != active_edge]
                
                temp_edge = Edge(-1, active_vertex, Vertex(-1, active_vertex.coords + array([.01,0])))
                active_edge = find_most_CCW(active_vertex, temp_edge, next_edge_candidates)
                active_edge.boundary = True
                #print "active edge: ", active_edge
                
                active_vertex = active_edge.other_vertex(active_vertex)
                active_vertex.bounadry = True
                #print "active vertex: ", active_vertex
                
                boundary_finder_i += 1
                
                #pdb.set_trace()
                
                continue
                
            next_edge_candidates = [edge for edge in self.Edge_list if edge.is_a_vertex(active_vertex) >=0 and edge != active_edge]
            active_edge = find_most_CW(active_vertex, active_edge, next_edge_candidates)
            active_edge.boundary = True
            #print "active edge: ", active_edge
            
            active_vertex = active_edge.other_vertex(active_vertex)
            active_vertex.bounadry = True
            #print "active vertex: ", active_vertex
            
            boundary_finder_i += 1
            
            #pdb.set_trace()
            if active_vertex ==find_most_SW(self.Vertex_list):
                return


if __name__ == '__main__':
    test_case = 5
    save_mesh = True
    save_mesh_filename = 'testmeshyue.pickle'
    
    if test_case == 1:
        vertex_list = [(0,0), (1,0), (2,0), (0,1), (1,1), (2,1), (1,2), (2,2)]
        edge_list = [(0,1), (1,2), (0,3), (1,4), (2,5), (3,4), (4,5), (4,6), (5,7), (6,7)]
        
    elif test_case == 2:
        vertex_list = [(18, 6), (0,0), (7,0), (13,0), (18,1), (3,5), (5,5), (11,8)]
        edge_list = [(1,2), (2,3), (3,4), (1,5), (2,6), (3,7), (0,4), (5,6), (6,7), (0,7)]
    
    elif test_case == 3:
        vertex_list = [(0,0), (4,0), (9,0), (13,0), (5,3), (8, 3), (3, 5), (10, 5)]
        edge_list = [(0,6), (0,1), (1,2), (2,3), (1,4), (2,5), (4,5), (6,7), (3,7), (4, 6), (5, 7)]
    
    elif test_case == 4:
        vertex_list = [(0,0), (4,0), (9,0), (13,0), (5,3), (8, 3), (3, 5), (10, 5), (0,8), (3,8), \
                    (10,8), (13,8), (-3,8), (-3,0), (16,8), (16,0), (-3,-2), (0,-2), (4,-2), (9,-2), \
                    (13,-2), (16,-2)]
                    
        edge_list = [(0,6), (0,1), (1,2), (2,3), (1,4), (2,5), (4,5), (6,7), (3,7), (4, 6),\
                (5, 7), (0,8), (6,9), (7,10), (3,11), (8,9), (9,10), (10,11), (8,12), (11,14),\
                (12,13), (14, 15), (13,16), (0,17), (18,1), (19,2), (20, 3), (15,21), (16,17), (17,18),\
                (18,19), (19,20), (20,21), (0,13), (3,15)]
                
    elif test_case == 5:
        dxf = DXFtoSegments.DXFGeometry('C:\cygwin64\home\Kerry\DXFtoMesh\Ampoule2.dxf')
        vertex_list,edge_list,bulge_list = dxf.cats2d_convert(len_scale=6)
        
    mesh = C2DMesh(vertex_list, edge_list, createC2Dmesh = False,  plotting = False, plotting_pause = 0.5)
    
    if save_mesh:
        with open(save_mesh_filename, 'w') as ff:
            pickle.dump(mesh, ff)


import DXFtoSegments
import pdb
import dxfgrabber
import pickle
import sys
from c2d_premesh_v5 import *
import os

'''Jeff's notes:
Changed "Nodes" to "nodes" for Pexpect
Added path control for working directory
'''

if 'darwin' in sys.platform:
    my_os = 'osx'
    import pexpect
elif 'linux' in sys.platform:
    my_os = 'linux'
    import pexpect
elif 'win32' in sys.platform:
    my_os = 'windows'
    import pexpect
elif 'cygwin' in sys.platform:
    my_os = 'cygwin'
    import pexpect
else:
    my_os = 'unknown:' + sys.platform
    import pexpect

    
def make_c2d_mesh(mesh, cats2d_path, working_dir='.'):
    # Get the current working directory (just in case '.' doesn't work)
    if working_dir == '.':
        working_dir = os.getcwd()

    mesh_filename = 'flow.mshc'
    meshplot_filename = 'mesh_plot.eps'
    
    if os.path.isfile(os.path.join(working_dir, 'flow.mshc')):
        overwrite_flowmshc = raw_input('A flow.mshc has been detected in {}.  Overwrite? (y/n): '.format(working_dir))
        
        if overwrite_flowmshc == 'y':
            os.remove(os.path.join(working_dir, 'flow.mshc'))
            
            if os.path.isfile(os.path.join(working_dir, 'mesh_plot.eps')):
                os.remove(os.path.join(working_dir, 'mesh_plot.eps'))
            
        else:
            mesh_filename = raw_input('Input new name for mesh file: ')
            
    
    directions = {'SOUTH':'s', 'NORTH':'n', 'EAST':'e', 'WEST':'w', 'South West': 'w', 'North West':'nw', 'South East': 'se', 'North East':'ne'}
    # Begins Cats2D process
    
    current_dir = os.getcwd()
    os.chdir(working_dir)
    print 'working directory: {}'.format(working_dir)
    cats2d = pexpect.spawn(cats2d_path, timeout=3, maxread=4000)
    fout = open('pexpect_log.txt', 'w')
    cats2d.logfile = fout
    #ssh_handle.logfile_read = sys.stdout
    #File 'flow.mshc' does not exist.
    #Try a new name? [y/N]: 
    i = cats2d.expect('name', timeout=3)
    b = cats2d.before
    a = cats2d.after
    cats2d.sendline('N')
    #Begin session with mesh generator? [Y/n]: Y
    i = cats2d.expect('[Y/n]')
    cats2d.sendline('Y')
    #Enter the Number of Regions, or 0 for array >= 0:
    i = cats2d.expect('array >= 0')
    cats2d.sendline(str(len(mesh.Region_list)))
    
    
    ## This is the part where it asks relations between one region to another.
    loop_i = 0
    print "Entering region relationship loop"
    while True:
        
        try:
            i = cats2d.expect('Enter the Neighbor Id', timeout = 1)
            
            #print cats2d
            
            
            #print "Start matching Neighbors!"
            
            b = cats2d.before.encode('utf8')
            a = cats2d.after.encode('utf8')
            b = b.translate(None, '\n\r?').strip(': ')
            neighbor_strings = b.split(' ')
            
            #card_dir = directions[neighbor_strings[-5]]
            region_id = int(neighbor_strings[-1])
            
            if 'NORTH' in neighbor_strings:
                card_dir = 'n'
            elif 'EAST' in neighbor_strings:
                card_dir = 'e'
            elif 'SOUTH' in neighbor_strings:
                card_dir = 's'
            elif 'WEST' in neighbor_strings:
                card_dir = 'w'
                
            #print "Made it here"
        
            neighbor_region = mesh.Region_list[region_id-1].Neighbors(card_dir)
            
            if neighbor_region != 0:
                #print "Made it to the if"
                neighbor_side = neighbor_region[1]
                cats2d.sendline(str(neighbor_region[0]))
                i = cats2d.expect('side', timeout = 3)
                #print "Which side?"
                cats2d.sendline(neighbor_side)
                print "Region " + str(region_id) + " has " + card_dir + " neighbor: Region " + str(neighbor_region[0]) + " (" + neighbor_side + ")"
            else:
                cats2d.sendline(str(neighbor_region))
                print "Region " + str(region_id) + " has no neighbor in the " + card_dir + " direction"
                    
            
        except pexpect.TIMEOUT:
            print "Jumping out of the loop!"
            break
    
    print "\n Entering vertex location loop"
    vertex_loop = 1
    while True:
        try:
            
            if vertex_loop > 1:
                cats2d.expect('value', timeout = 3)
            
            b = cats2d.before.encode('utf8')
            #pdb.set_trace()
            
            if b.find('Region')>0:            
                vertex_strings = b.translate(None, '\n\r?').split("="*79)
                region_id = int(filter(str.isdigit, vertex_strings[-2]))
                
                #pdb.set_trace()
                
            corner_string = b.split("***")[-2]
            
            if corner_string.rfind('South West') > 0:
                corner_dir = 'sw'
            elif corner_string.rfind('North West') > 0:
                corner_dir = 'nw'
            elif corner_string.rfind('South East') > 0:
                corner_dir = 'se'
            elif corner_string.rfind('North East') > 0:
                corner_dir = 'ne'
            
            #if region_id >=10:
            #    pdb.set_trace()
                
                        
            x_corner = str(mesh.Region_list[region_id - 1].Vertices(corner_dir).coords[0])
            y_corner = str(mesh.Region_list[region_id - 1].Vertices(corner_dir).coords[1])
            
            cats2d.sendline(x_corner)
            cats2d.expect('y value', timeout = 3)
            cats2d.sendline(y_corner)
            
            print "The " + corner_dir + " corner of Region " + str(region_id) + " is at (" + x_corner + "," + y_corner + ")"
            vertex_loop += 1
            #print vertex_loop
            
        except pexpect.TIMEOUT:
            print "Jumping out of the loop!"
            break
    
    print "\nEntering discretization loop"
    discretization_loop = 1
    while True:
        try:
            if discretization_loop > 1:
                cats2d.expect('elements', timeout = 3)
            
            b = cats2d.before.encode('utf8').translate(None, '\n\r?')
            bc = b.split("="*79)
            
            #pdb.set_trace()
            
            if b.find('Enter the Number of')>0:            
                cats2d.sendline('1')
                cats2d.expect('Enter the Number of', timeout = 3)
                cats2d.sendline('1')
                
            discretization_loop += 1
            print discretization_loop
            
        except pexpect.TIMEOUT:
            print "Jumping out of the loop!"
            break
    
    #pdb.set_trace()
    
    cats2d.expect('M E S H   D R I V E R')
    cats2d.sendline('Quit')
    cats2d.expect('Save changes to file?')
    cats2d.sendline('Y')
    cats2d.expect('Enter the restart file name:')
    cats2d.sendline(mesh_filename)
    print 'Mesh file written as ' + mesh_filename
    cats2d.expect('C A T S   2 D')
    cats2d.sendline('Post Process')
    cats2d.expect('P O S T   P R O C E S S O R')
    cats2d.sendline('MeshPlot')
    cats2d.expect('Enter the plot file name:')
    cats2d.sendline(meshplot_filename)
    print 'Mesh plot written as ' + meshplot_filename
    cats2d.expect('Plot nodes?')
    cats2d.sendline('N')
    cats2d.expect('P O S T   P R O C E S S O R')
    cats2d.sendline('Return')
    cats2d.expect('C A T S   2 D')
    cats2d.sendline('Quit')
    os.chdir(current_dir)

    
    fout.close() 
    print 'Complete'
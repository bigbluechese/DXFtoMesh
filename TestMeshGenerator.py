import DXFtoSegments
import c2d_premesh_v5
import pexpect_c2dmesh_v2

dxf = DXFtoSegments.DXFGeometry('./DXFTests/Ampoule2.dxf')
vertex_list,edge_list,bulge_list = dxf.cats2d_convert(len_scale=6)

mesh = c2d_premesh_v5.C2DMesh(vertex_list, edge_list, plotting = False, plotting_pause = 0.1)
cats2d_path = 'cats2d.x'
pexpect_c2dmesh_v2.make_c2d_mesh(mesh, cats2d_path, working_dir=dxf.work_dir)
import pytools as pt
import numpy as np
import sys

run = sys.argv[1]

#path_bulk = '/wrk/users/dubart/300_test/proc_test/'+run+'/'
#path_bulk = '/wrk/users/dubart/300_test/'+run+'/bulk/'
#path_bulk = '/wrk-vakka/users/dubart/diff_test/dmumu/'+run+'/bulk/'
#path_bulk = '/wrk/users/dubart/diff_test/dmumu/diff/'+run+'/CFL1.0/bulk/'
path_bulk = '/wrk-vakka/users/markusb/weirddiffusion/Losscone_may30/'
path_save = '/wrk-vakka/users/dubart/analysis/subgrid/300_test/Green_diff/data/mu/'+run+'/'

RE = 6371e3 # m

# Source data file
bulkname = "bulk.0000010.vlsv"
    
#Open file
f = pt.vlsvfile.VlsvReader(path_bulk+bulkname)

xmin = f.read_parameter('xmin')
ymin = f.read_parameter('ymin')
zmin = f.read_parameter('zmin')

xmax = f.read_parameter('xmax')
ymax = f.read_parameter('ymax')
zmax = f.read_parameter('zmax')

coords  = [(xmin+xmax)/2.0,(ymin+ymax)/2.0,(zmin+zmax)/2.0]

#CellID = f.read_variable("CellID")
CellID = int(f.get_cellid(coords))
print(CellID)

cell_candidates = f.read(mesh = "SpatialGrid", tag = "CELLSWITHBLOCKS")
# Read in the coordinates of the cells:
cell_candidate_coordinates = [f.get_cell_coordinates(cell_candidate) for cell_candidate in cell_candidates]

# Read in the cell's coordinates:
pick_cell_coordinates = f.get_cell_coordinates(CellID)
if len(cell_candidates) == 0:
    print("No velocity distribution data found in this file!")
    sys.exit()

# Find the nearest:
from operator import itemgetter
norms = np.sum((cell_candidate_coordinates - pick_cell_coordinates)**2, axis=-1)**(1./2)
norm, i = min((norm, idx) for (idx, norm) in enumerate(norms))

# Get the cell id:
cid = cell_candidates[i]
print("PLOTTED CELL ID: " + str(cid))
coords = f.get_cell_coordinates(cid)
print(str(coords/RE))

# check if velocity space exists in this cell
if f.check_variable('fSaved'): #restart files will not have this value                                                                      
    if f.read_variable('fSaved',cid) != 1.0:
        print('No vdf in this cell')

print("CHECKED")

np.save(path_save+'CellID_'+run+'.npy',cid)
print('Saved '+path_save+'CellID_'+run+'.npy')

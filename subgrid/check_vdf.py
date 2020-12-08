import pytools as pt
import numpy as np
import sys

run = sys.argv[1]

path_bulk = '/wrk/users/dubart/'+run+'km/maxwellian_periodic/'
#path_bulk = '/wrk/users/dubart/'+run+'km/Shock/'
path_save = '/wrk/users/dubart/'+run+'km/data/'

RE = 6371e3 # m

if len(sys.argv)==4: # Starting and end bulk given
    timetot = range(int(sys.argv[2]), int(sys.argv[3]), 1)
else: # Only starting bulk given, generate one bulk
    timetot = range(int(sys.argv[2]), int(sys.argv[2])+1, 1)
for j in timetot:
    # Source data file
    bulkname = "bulk."+str(j).rjust(7,'0')+".vlsv"
    print(bulkname)
    
    #Open file
    #f = pt.vlsvfile.VlsvReader(path_VDF+distribname)
    f = pt.vlsvfile.VlsvReader(path_bulk+bulkname)

    xmin = f.read_parameter('xmin')
    ymin = f.read_parameter('ymin')
    zmin = f.read_parameter('zmin')

    xmax = f.read_parameter('xmax')
    ymax = f.read_parameter('ymax')
    zmax = f.read_parameter('zmax')
    
    coord_mid = [(xmin+xmax)/2.0,(ymin+ymax/2.0),(zmin+zmax)/2.0]

    CellID = int(f.get_cellid(coord_mid))
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
    print(coords)
    np.save(path_save+'CellID_box.npy',cid)

    # check if velocity space exists in this cell
    if f.check_variable('fSaved'): #restart files will not have this value                                                                      
        if f.read_variable('fSaved',cid) != 1.0:
            print('No vdf in this cell')

    print("CHECKED")



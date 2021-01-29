import pytools as pt
import numpy as np
import sys
from multiprocessing import Pool

numproc = 20

run = sys.argv[1]

bulkStart = int(sys.argv[2])
bulkEnd   = int(sys.argv[3])

#path_bulk = '/wrk/users/dubart/'+run+'km/maxwellian_periodic/'
path_bulk = '/wrk/users/dubart/'+run+'km/diffusion/bulk/'
path_save = '/wrk/users/dubart/analysis/subgrid/300_test/diff_test/data/'

RE = 6371e3 # m

CellID = int(np.load(path_save+'CellID_'+run+'.npy'))

def get_variables(step):

    bulkname = 'bulk.'+str(step).rjust(7,'0')+'.vlsv'
    f = pt.vlsvfile.VlsvReader(path_bulk+bulkname)

    beta      = f.read_variable('vg_beta',CellID)
    beta_par  = f.read_variable('vg_beta_parallel',CellID)
    beta_perp = f.read_variable('vg_beta_perpendicular',CellID)

    Temperature = f.read_variable('vg_temperature',CellID)
    T_par       = f.read_variable('vg_t_parallel',CellID)
    T_perp      = f.read_variable('vg_t_perpendicular',CellID)
    T_aniso     = f.read_variable('vg_t_anisotropy',CellID)

    return [beta,beta_par,beta_perp,Temperature,T_par,T_perp,T_aniso]


pool = Pool(numproc)
data = pool.map(get_variables,range(bulkStart,bulkEnd+1))
pool.terminate()

data = np.array(data)

beta      = data[:,0]
beta_par  = data[:,1]
beta_perp = data[:,2]

Temperature = data[:,3]
T_par       = data[:,4]
T_perp      = data[:,5]
T_aniso     = data[:,6]

if path_bulk == '/wrk/users/dubart/900km/diffusion/bulk/':
    run = '900_diff'

np.save(path_save+'/beta_'+run+'.npy',beta)
np.save(path_save+'/beta_par_'+run+'.npy',beta_par) 
np.save(path_save+'/beta_perp_'+run+'.npy',beta_perp)

np.save(path_save+'/Temperature_'+run+'.npy',Temperature)
np.save(path_save+'/T_par_'+run+'.npy',T_par)    
np.save(path_save+'/T_perp_'+run+'.npy',T_perp)
np.save(path_save+'/T_aniso_'+run+'.npy',T_aniso)

print('Saved')

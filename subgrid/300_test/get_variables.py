import pytools as pt
import numpy as np
import sys
from multiprocessing import Pool

numproc = 20

run = sys.argv[1]

bulkStart = 0000
bulkEnd   = 1000

path_bulk = '/wrk/users/dubart/300_test/'+run+'/bulk/'
path_save = '/wrk/users/dubart/analysis/subgrid/300_test/data/'

RE = 6371e3 # m

loc = ['mid','right','left','top_mid','bot_mid','top_r','bot_r','top_l','bot_l']

CellIDs = np.zeros(len(loc))

for i in range(0,len(loc)):

    CellIDs[i] = np.load(path_save+'CellID_'+loc[i]+'.npy')


def get_variables(step):

    bulkname = 'bulk.'+str(step).rjust(7,'0')+'.vlsv'
    f = pt.vlsvfile.VlsvReader(path_bulk+bulkname)

    beta      = f.read_variable('vg_beta',CellIDs)
    beta_par  = f.read_variable('vg_beta_parallel',CellIDs)
    beta_perp = f.read_variable('vg_beta_perpendicular',CellIDs)

    Temperature = f.read_variable('vg_temperature',CellIDs)
    T_par       = f.read_variable('vg_t_parallel',CellIDs)
    T_perp      = f.read_variable('vg_t_perpendicular',CellIDs)
    T_aniso     = f.read_variable('vg_t_anisotropy',CellIDs)

    return [beta,beta_par,beta_perp,Temperature,T_par,T_perp,T_aniso]


pool = Pool(numproc)
data = pool.map(get_variables,range(bulkStart,bulkEnd+1))
pool.terminate()

data = np.array(data)

beta      = data[:,0,:]
beta_par  = data[:,1,:]
beta_perp = data[:,2,:]

Temperature = data[:,3,:]
T_par       = data[:,4,:]
T_perp      = data[:,5,:]
T_aniso     = data[:,6,:]

np.save(path_save+run+'/beta.npy',beta)
np.save(path_save+run+'/beta_par.npy',beta_par) 
np.save(path_save+run+'/beta_perp.npy',beta_perp)

np.save(path_save+run+'/Temperature.npy',Temperature)
np.save(path_save+run+'/T_par.npy',T_par)    
np.save(path_save+run+'/T_perp.npy',T_perp)
np.save(path_save+run+'/T_aniso.npy',T_aniso)

print('Saved')

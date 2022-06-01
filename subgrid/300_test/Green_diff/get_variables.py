import pytools as pt
import numpy as np
import sys
from multiprocessing import Pool

numproc = 20

A   = sys.argv[1]
run = sys.argv[2]

bulkStart = int(sys.argv[3])
bulkEnd   = int(sys.argv[4])

Dmumu = sys.argv[5]
CFL   = sys.argv[6]

path_bulk = '/wrk/users/dubart/diff_test/proc_test/'+run+'/'
#path_bulk = '/wrk/users/dubart/300_test/900km/diffusion/restart/bulk/'
path_save = '/wrk/users/dubart/analysis/subgrid/300_test/Green_diff/data/'

RE = 6371e3 # m

#CellID = int(np.load(path_save+'CellID_diffusion.npy'))
CellID = 1

def get_variables(step):

    if step == 0:
        bulkname = 'initial-grid.0000000.vlsv'
    else:
        bulkname = 'bulk.'+str(step).rjust(7,'0')+'.vlsv'
   
    f = pt.vlsvfile.VlsvReader(path_bulk+bulkname)

    rho         = f.read_variable('vg_rho',CellID)
    Temperature = f.read_variable('vg_temperature',CellID)
    Tpara       = f.read_variable('vg_t_parallel',CellID)
    Tperp       = f.read_variable('vg_t_perpendicular',CellID)
    beta        = f.read_variable('vg_beta',CellID)
    beta_para   = f.read_variable('vg_beta_parallel',CellID)
    beta_perp   = f.read_variable('vg_beta_perpendicular',CellID)

    return [rho,Temperature,Tpara,Tperp,beta,beta_para,beta_perp]


pool = Pool(numproc)
data = pool.map(get_variables,range(bulkStart,bulkEnd+1))
pool.terminate()

data = np.array(data)

rho         = data[:,0]
Temperature = data[:,1]
Tpara       = data[:,2]
Tperp       = data[:,3]
beta        = data[:,4]
beta_para   = data[:,5]
beta_perp   = data[:,6]


np.save(path_save+'A_'+A+'_rho_'+run+'_D'+Dmumu+'_CFL'+CFL+'.npy',rho)
np.save(path_save+'A_'+A+'_T_'+run+'_D'+Dmumu+'_CFL'+CFL+'.npy',Temperature)
np.save(path_save+'A_'+A+'_Tpara_'+run+'_D'+Dmumu+'_CFL'+CFL+'.npy',Tpara)
np.save(path_save+'A_'+A+'_Tperp_'+run+'_D'+Dmumu+'_CFL'+CFL+'.npy',Tperp)
np.save(path_save+'A_'+A+'_beta_'+run+'_D'+Dmumu+'_CFL'+CFL+'.npy',beta)
np.save(path_save+'A_'+A+'_beta_para_'+run+'_D'+Dmumu+'_CFL'+CFL+'.npy',beta_para)
np.save(path_save+'A_'+A+'_beta_perp_'+run+'_D'+Dmumu+'_CFL'+CFL+'.npy',beta_perp)

print(path_save+'A_'+A+'_*_'+run+'_D'+Dmumu+'_CFL'+CFL+'.npy')

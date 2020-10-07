import numpy as np
import os,sys
from multiprocessing import Pool

numproc=40

path_HYDROS = '/wrk/users/dubart/analysis/hydros/HYDROS/'
path_output = '/wrk/users/dubart/analysis/hydros/Lucile/output/'
path_file   = '/wrk/users/dubart/analysis/hydros/Lucile/write/'

theta   = np.arange(70.0,90.1,0.1)
w_start = 1.0

def run_hydros(step):

    filename = "Lucile_"+str(w_start)+'_'+str(round(theta[step],1))

    os.system("python "+path_HYDROS+"hydros.py "+path_file+filename)

    return

pool = Pool(numproc)
pool.map(run_hydros,range(0,len(theta)))

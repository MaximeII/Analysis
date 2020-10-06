import numpy as np
import os,sys

run = sys.argv[1]

path_HYDROS = '/wrk/users/dubart/analysis/hydros/HYDROS/'
path_file   = '/wrk/users/dubart/analysis/hydros/'+run+'/'
path_output = '/wrk/users/dubart/analysis/hydros/output/'+run+'/'

k  = np.linspace(0.01,1.2,121)

for i in range(0,len(k)):

    os.system("python "+path_HYDROS+"hydros.py "+path_file+run+'_beta2_'+str(round(k[i],2))+'_45')


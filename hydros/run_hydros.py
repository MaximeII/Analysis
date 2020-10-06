import numpy as np
import os,sys

run = sys.argv[1]

path_HYDROS = '/wrk/users/dubart/analysis/hydros/HYDROS/'
path_file   = '/wrk/users/dubart/analysis/hydros/'+run+'/'
path_output = '/wrk/users/dubart/analysis/hydros/output/'+run+'/'

kpar  = np.linspace(0.01,1.5,151)

if run == 'BCG' or run == 'BGA':

    os.system("python "+path_HYDROS+"hydros.py "+path_file+run+'_krange_beta2')

else: 

    for i in range(0,len(kpar)):

        os.system("python "+path_HYDROS+"hydros.py "+path_file+run+'_beta2_'+str(round(kpar[i],2)))

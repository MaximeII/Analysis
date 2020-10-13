import numpy as np
import sys
from multiprocessing import Pool

numproc = 40

path_HYDROS = '/wrk/users/dubart/analysis/hydros/HYDROS/'
path_output = '/wrk/users/dubart/analysis/hydros/Lucile/output/'
path_file   = '/wrk/users/dubart/analysis/hydros/Lucile/write/'


wavemode = '2' #1: use kperp, kpar; =2: use k, theta

beta_para = 2.5
Taniso    = 1.0

#Constants (SI units)
q        = 1.60217e-19
mu0      = 1.256637061e-6
c        = 299792458
mp       = 1.6726219e-27
kB       = 1.38064852e-23
me       = 9.1093826e-31
epsilon0 = 8.85418781762e-12

tau = 1.0

w_start = 0.01

scan_type = '1d'

theta = np.arange(0.0,40.1,0.1)

def write_hydros(step):

    filename = "Lucile_"+str(w_start)+'_'+str(round(theta[step],1))
    
    template   = open(path_HYDROS+'1dscan', 'r')
    HYDROSfile = open(path_file+filename, 'w')
    for line in template:
        if line[:10] == "wavevector":
            HYDROSfile.write('wavevector_mode='+wavemode+'\n')
        elif line[:2] == "k=":
            HYDROSfile.write('k= 0.01\n')
        elif line[:6] == "theta=":
            HYDROSfile.write('theta= '+str(round(theta[step],1))+'\n')
        elif line[:5] == "beta=":
            HYDROSfile.write('beta= '+str(beta_para)+'\n')
        elif line[:4] == "tau=":
            HYDROSfile.write('tau= '+str(tau)+'\n')
        elif line[:11] == "Tpar_Tperp=":
            HYDROSfile.write('Tpar_Tperp= '+str(round(1/Taniso,2))+'\n')
        elif line[:6] == "start=":
            HYDROSfile.write('start= ['+str(w_start)+' , -0.]\n')
        elif line[:14] == "normalization=":
            HYDROSfile.write('normalization=1\n')
        elif line[:10] == "scan_type=":
            HYDROSfile.write("scan_type= "+scan_type+"\n")
        elif line[:8] == "scanvar=":
            HYDROSfile.write("scanvar= 'k'\n") 
        elif line[:11] == "scan_range=":
            HYDROSfile.write('scan_range= (0.01,5.0)\n')
        elif line[:5] == "nval=":
            HYDROSfile.write('nval= 500\n')
        elif line[:9] == "scanvar2=":
            HYDROSfile.write("scanvar2= 'theta'\n")    
        elif line[:12] == "scan_range2=":
            HYDROSfile.write('scan_range2= (88.0,90.0)\n')
        elif line[:6] == "nval2=":
            HYDROSfile.write('nval2= 5\n')
        elif line[:11] == "resultfile=":
            HYDROSfile.write("resultfile= '"+path_output+filename+".log'\n")
        else:
            HYDROSfile.write(line)
    HYDROSfile.close() 
    print("Wrote "+path_file+filename)
   
    return

pool = Pool(numproc)
print('Writing ...')
pool.map(write_hydros,range(0,len(theta)))

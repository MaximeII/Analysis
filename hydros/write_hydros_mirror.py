import numpy as np
import sys

run = sys.argv[1]

path_HYDROS = '/wrk/users/dubart/analysis/hydros/HYDROS/'
path_file   = '/wrk/users/dubart/analysis/hydros/'+run+'/'
path_output = '/wrk/users/dubart/analysis/hydros/output/'+run+'/'

wavemode = '2' #1: use kperp, kpar; =2: use k, thetaa

k     = np.linspace(0.01,1.2,121)
theta = 45

if run == 'BCQ':
    Taniso = 2.0390744819641364
    Tpara  = 5564061.811271753     # K
    Tperp  = 11291411.958148913    # K
    T      = 9382295.242523193     # K
    B      = 1.951962364912315e-08 # T
    Beta   = 2.5838120374104614
    n      = 2934092.7131413473    # m-3
    r_L    = 230913.6169381046     # m
    d_i    = 132937.78719913246    # m
    dx     = 300e3
elif run == 'BGA':
    Taniso = 2.684076196248551
    Tpara  = 4695840.183918735      # K
    Tperp  = 12431821.920983285     # K
    T      = 9853161.34196177       # K
    B      = 1.9874046310499582e-08 # T
    Beta   = 2.7173036160675155
    n      = 2661472.5324297952     # m-3
    r_L    = 237973.1450922526      # m
    d_i    = 139580.37722259495     # m
    dx     = 600e3
elif run == 'BCG':
    Taniso = 3.334012229418161
    Tpara  = 4320021.260923261      # K
    Tperp  = 13273600.20811248      # K
    T      = 10289073.892382741     # K
    B      = 1.6404192010995193e-08 # T
    Beta   = 4.182207359983115
    n      = 2850283.4332915875     # m-3
    r_L    = 297910.8906913871      # m
    d_i    = 134878.06835706628     # m
    dx     = 900e3

#Constants (SI units)
q        = 1.60217e-19
mu0      = 1.256637061e-6
c        = 299792458
mp       = 1.6726219e-27
kB       = 1.38064852e-23
me       = 9.1093826e-31
epsilon0 = 8.85418781762e-12

beta_para = n * kB * Tpara / ( B**2 / (2*mu0)) *2
#tau = 1.0
tau = Tpara * 8.621738e-5 / 1.0
v_A  = B / np.sqrt(mu0 * mp * n)
w_ci = q * B / mp

for i in range(0,len(k)):

    template   = open(path_HYDROS+'1dscan', 'r')
    HYDROSfile = open(path_file+run+'_beta2_'+str(round(k[i],2))+'_'+str(theta), 'w')
    for line in template:
        if line[:10] == "wavevector":
            HYDROSfile.write('wavevector_mode='+wavemode+'\n')
        elif line[:3] == "k= ":
            HYDROSfile.write('k= '+str(round(k[i],2))+'\n')
        elif line[:6] == "theta=":
            HYDROSfile.write('theta= '+str(theta)+'\n')
        elif line[:5] == "beta=":
            HYDROSfile.write('beta= '+str(round(beta_para,2))+'\n')
        elif line[:4] == "tau=":
            HYDROSfile.write('tau= '+str(round(tau,2))+'\n')
        elif line[:11] == "Tpar_Tperp=":
            HYDROSfile.write('Tpar_Tperp= '+str(round(1/Taniso,2))+'\n')
        elif line[:6] == "start=":
            HYDROSfile.write('start= [0.0 , -0.]\n')
        elif line[:14] == "normalization=":
            HYDROSfile.write('normalization=1\n')
        elif line[:8] == "scanvar=":
            HYDROSfile.write("scanvar= 'k'\n") 
        elif line[:11] == "scan_range=":
            HYDROSfile.write('scan_range= ('+str(round(k[i],2))+','+str(round(k[i],2))+')\n')
        elif line[:5] == "nval=":
            HYDROSfile.write('nval= 1\n')
        elif line[:11] == "resultfile=":
            HYDROSfile.write("resultfile= '"+path_output+run+"_beta2_"+str(round(k[i],2))+'_'+str(theta)+".log'\n")    
        else:
            HYDROSfile.write(line)
    HYDROSfile.close() 
    print('Wrote '+run+'_beta2_'+str(round(k[i],2))+'_'+str(theta))  

import numpy as np
import scipy.integrate
import matplotlib
import matplotlib.pyplot as plt
import sys
import pytools as pt

path_fig   = '/wrk-vakka/users/dubart/analysis/subgrid/300_test/Green_diff/fig/mu/proc_test/'
path_files = '/wrk-vakka/users/dubart/diff_test/proc_test/mu_files/'
path_bulk  = '/wrk-vakka/users/dubart/diff_test/proc_test/bulk/'

if len(sys.argv) == 3: # Starting and end frames given
    timetot = range(int(sys.argv[1]), int(sys.argv[2]),1)
else: # Only starting frame given, generate one frame
    timetot = range(int(sys.argv[1]), int(sys.argv[1])+1, 1)
for step in timetot:

    fig, axs = plt.subplots(2,1,sharex=True)

    if step == 0:
        bulkname = "initial-grid.0000000.vlsv"
    else:
        bulkname = "bulk."+str(step).rjust(7,'0')+".vlsv"

    cid = 1

    # Get/plot fmu data
    data_fmu = np.loadtxt(path_files+'muv_array_'+str(step).rjust(7,'0')+'.txt')
    
    murange = np.linspace(-1.0, 1.0,data_fmu.shape[1])
    dmu     = abs(murange[0] - murange[1])
    fmumin  = 0.0
    fmumax  = 0.4
    fmu = np.mean(data_fmu,axis=0)    

    axs[0].plot(murange,fmu,linewidth=2,color='black')
    axs[0].set_xlim(-1.0,1.0)
    axs[0].set_ylim(fmumin,fmumax)
    axs[0].set_ylabel('f($\\mu$)',fontsize=20,weight='bold')
    
    # Compute dfdmu and dfdmumu
    dfdmu   = np.empty(len(fmu))
    dfdmumu = np.empty(len(fmu))
        # dfdmu
    dfdmu[0]  = (fmu[1] - fmu[0]) / dmu         
    dfdmu[1:-1] = (fmu[2:] - fmu[:-2]) / (2*dmu)  
    dfdmu[-1] = (fmu[-1] - fmu[-2]) / dmu
        # dfdmumu
    dfdmumu[0]    = (fmu[2] - 2*fmu[1] + fmu[0]) / dmu**2
    dfdmumu[1:-1] = (fmu[2:] - 2*fmu[1:-1] + fmu[:-2]) / dmu**2
    dfdmumu[-1]   = (fmu[-1] - 2*fmu[-2] + fmu[-3]) / dmu**2

    # Get dfdt_mu
    data_dfdt = np.loadtxt(path_files+'dfdt_mu_array_'+str(step).rjust(7,'0')+'.txt')
    dfdt_mu   = np.mean(data_dfdt,axis=0)

    Dmumu = np.empty(len(fmu))

    # Compute Dmumu
    Matrix = np.diag(dfdmumu,0) + np.diag(dfdmu[:-1]/(2*dmu),1) + np.diag(-dfdmu[1:]/(2*dmu),-1)
    try:
       InvM   = np.linalg.inv(Matrix)
       Dmumu  = np.matmul(InvM,dfdt_mu)
    except Exception as e:
       print('step = '+str(step)+': '+str(e)) 

   
    axs[1].plot(murange,Dmumu,linewidth=2,color='black')
    axs[1].set_xlim(-1.0,1.0)
    axs[1].set_ylim(0.0,0.1)
    axs[1].set_xlabel('$\\mu$',fontsize=20)
    axs[1].set_ylabel('D$_{\\mu\\mu}$',fontsize=20,weight='bold')

    plt.savefig(path_fig+'Dmumu_proc_test_'+str(step).rjust(7,'0')+'.png',dpi=300)
    print(path_fig+'Dmumu_proc_test_'+str(step).rjust(7,'0')+'.png')
    plt.close()

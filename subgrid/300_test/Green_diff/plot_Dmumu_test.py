import pytools as pt
import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt

run=sys.argv[1]

path_save = '/wrk/users/dubart/analysis/subgrid/300_test/Green_diff/data/mu/'+run+'/'
path_bulk = '/wrk/users/dubart/diff_test/dmumu/'+run+'/bulk/'
path_fig  = '/wrk-vakka/users/dubart/analysis/subgrid/300_test/Green_diff/fig/mu/'

murange  = np.linspace(-1.0,1.0,101)
time = np.arange(0.0,12.0,0.01)

fmu   = np.load(path_save+'fmu_1D_test.npy')
Dmumu = np.load(path_save+'Dmumu_1D_test.npy')

for i in range(0,len(time)):

    fig, axs = plt.subplots(2,1,sharex=True)

    axs[0].set_xlim(-1,1)
    axs[1].set_xlim(-1,1)

    axs[0].set_ylim(0.0,0.3)
    axs[1].set_ylim(0.0,0.02)

    axs[0].set_ylabel('f($\\mu$)',fontsize=20,weight='bold')

    axs[1].set_xlabel('$\\mu$',fontsize=20,weight='bold')
    axs[1].set_ylabel('D$_{\\mu\\mu}$ [s$^{-1}$]',fontsize=20,weight='bold')

    axs[0].plot(murange,fmu[:,i],linewidth=2,color='black')
    axs[1].plot(murange,Dmumu[:,i],linewidth=2,color='black')

    plt.suptitle('t = '+str(round(time[i],3))+'s',fontsize=18,weight='bold')
    plt.savefig(path_fig+'Dmumu_1D_test_'+str(i).rjust(5,'0')+'.png',dpi=300)
    print(path_fig+'Dmumu_1D_test_'+str(i).rjust(5,'0')+'.png')
    plt.close()

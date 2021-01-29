import pytools as pt
import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt

numproc = 20

run = ['300','900','900_diff']

bulkStart = int(sys.argv[1])
bulkEnd   = int(sys.argv[2])

path_save = '/wrk/users/dubart/analysis/subgrid/300_test/diff_test/data/'
path_fig  = '/wrk/users/dubart/analysis/subgrid/300_test/diff_test/fig/'

RE = 6371e3 # m

fig, axs = plt.subplots(2,1,sharex=True)

for i in range(0,len(run)):

    path_bulk = '/wrk/users/dubart/'+run[i]+'km/maxwellian_periodic/'

    if i == 2:

        path_bulk = '/wrk/users/dubart/'+run[i]+'km/diffusion/bulk/'
    

    beta      = np.load(path_save+'/beta_'+run[i]+'.npy')
    beta_par  = np.load(path_save+'/beta_par_'+run[i]+'.npy') 
    beta_perp = np.load(path_save+'/beta_perp_'+run[i]+'.npy')
    
    Temperature = np.load(path_save+'/Temperature_'+run[i]+'.npy')
    T_par       = np.load(path_save+'/T_par_'+run[i]+'.npy')    
    T_perp      = np.load(path_save+'/T_perp_'+run[i]+'.npy')
    T_aniso     = np.load(path_save+'/T_aniso_'+run[i]+'.npy')

    time = np.linspace(bulkStart/2,bulkEnd/2,len(beta_par))

    axs[0].plot(time,beta_par,linewidth = 2,label = run[i])
    axs[1].plot(time,T_aniso,linewidth = 2)

axs[0].set_ylabel('$\\beta_{\\parallel}$',fontsize = 15)
axs[0].legend(bbox_to_anchor=(-0.1,0.95,1.2,0.2),loc='upper right',ncol = 3, mode = 'expand',fontsize = 7)

axs[1].set_ylabel('$T_{\\bot}/T_{\\parallel}$',fontsize = 15)
axs[1].set_xlabel('Time (s)')
axs[1].set_xlim(0.0,bulkEnd/2)

plt.savefig(path_fig+'Variables_diff.png',dpi=300)
print(path_fig+'Variables_diff.png')

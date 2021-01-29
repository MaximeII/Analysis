import pytools as pt
import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt

numproc = 20

run = sys.argv[1]

bulkStart = int(sys.argv[2])
bulkEnd   = int(sys.argv[3])

path_bulk = '/wrk/users/dubart/300_test/'+run+'/bulk/'
path_save = '/wrk/users/dubart/analysis/subgrid/300_test/data/'
path_fig  = '/wrk/users/dubart/analysis/subgrid/300_test/fig/'+run+'/'

RE = 6371e3 # m

loc = ['mid','right','left','top_mid','bot_mid','top_r','bot_r','top_l','bot_l']

beta      = np.load(path_save+run+'/beta.npy')
beta_par  = np.load(path_save+run+'/beta_par.npy') 
beta_perp = np.load(path_save+run+'/beta_perp.npy')

Temperature = np.load(path_save+run+'/Temperature.npy')
T_par       = np.load(path_save+run+'/T_par.npy')    
T_perp      = np.load(path_save+run+'/T_perp.npy')
T_aniso     = np.load(path_save+run+'/T_aniso.npy')

time = np.linspace(bulkStart,bulkEnd/2,beta.shape[0])

fig, axs = plt.subplots(2,1,sharex=True)

for i in range(0,len(loc)):

    axs[0].plot(time,beta_par[:,i],linewidth = 2,label = loc[i])
    axs[1].plot(time,T_aniso[:,i],linewidth = 2)

axs[0].set_ylabel('$\\beta_{\\parallel}$',fontsize = 15)
axs[0].legend(bbox_to_anchor=(-0.1,0.95,1.2,0.2),loc='upper right',ncol = 9, mode = 'expand',fontsize = 7)

axs[1].set_ylabel('$T_{\\bot}/T_{\\parallel}$',fontsize = 15)
axs[1].set_xlabel('Time (s)')
axs[1].set_xlim(0.0,bulkEnd/2)

fig.suptitle(run,fontsize = 15)

plt.savefig(path_fig+'Variables_'+run+'.png',dpi=300)
print(path_fig+'Variables_'+run+'.png')

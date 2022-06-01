import numpy as np
import matplotlib.pyplot as plt
import sys
import pytools as pt

path_fig  = '/wrk/users/dubart/analysis/subgrid/300_test/Green_diff/fig/evaluateA/'
path_data = '/wrk/users/dubart/analysis/subgrid/300_test/Green_diff/data/'
path_bulk = '/wrk/users/dubart/diff_test/proc_test/mu_bulk/'

f = pt.vlsvfile.VlsvReader(path_bulk+'initial-grid.0000000.vlsv')

Tpara = f.read_variable('vg_t_parallel')
Tperp = f.read_variable('vg_t_perpendicular')

Ainit = Tpara/Tperp

dt    = 0.01
Dmumu = [0.1,0.05,0.01,0.005,0.001]

CFL   = ['0.05','0.1','0.3','0.5','0.7']

Dleg   = ['0.1','0.05','0.01','0.005','0.001']
#var    = '$D_{\\mu\\mu}$'
var    = 'CFL'
panels = ['(a)','(b)','(c)','(d)','(e)']


fig, axs = plt.subplots(5,1,sharex = True)

tmax = 12.0
tmin = 0.0

Amin = Ainit
Amax = 1.0

axs[0].set_xlim(tmin,tmax)
axs[0].set_ylim(Amin,Amax)

axs[1].set_xlim(tmin,tmax)
axs[1].set_ylim(Amin,Amax)

axs[2].set_xlim(tmin,tmax)
axs[2].set_ylim(Amin,Amax)
axs[2].set_ylabel('$T_\\parallel/T_\\bot$',fontsize=15)

axs[3].set_xlim(tmin,tmax)
axs[3].set_ylim(Amin,Amax)

axs[4].set_xlim(tmin,tmax)
axs[4].set_ylim(Amin,Amax)
axs[4].set_xlabel('Time (s)',fontsize=15)

#for i in range(0,len(Dmumu)):
#    A = np.load(path_data+'A_'+"{:.2f}".format(Ainit)+'_D'+str(Dmumu[i])+'_dt0.001.npy')
#    trange = np.arange(tmin,tmax,0.001)    
#
#    axs[i].plot(trange,A,color='black',label='Theory')
#
#    try:
#
#        Tpara = np.load(path_data+'A_0.3_Tpara_mu_bulk_D'+str(Dmumu[i])+'_CFL0.1.npy')
#        Tperp = np.load(path_data+'A_0.3_Tperp_mu_bulk_D'+str(Dmumu[i])+'_CFL0.1.npy')
#           
#        Asim = Tpara/Tperp
#        tsim = np.arange(tmin,len(Asim)*dt,dt)
#
#        axs[i].plot(tsim,Asim,color='red',linestyle = '--',label='Simulation')
#   
#    except IOError:
#        print('No file found')
#        continue
#
#    if i == 0:
#        axs[i].legend(loc='best',bbox_to_anchor=(0.0,1.1,0.5,0.5),ncol=2)
#
#    axs[i].text(8.8,0.75,panels[i]+' '+var+' = '+Dleg[i],fontsize=12)


for i in range(0,len(CFL)):
    A = np.load(path_data+'A_'+"{:.2f}".format(Ainit)+'_D1.0_dt0.01.npy')
    trange = np.arange(tmin,tmax,0.01)    

    axs[i].plot(trange,A,color='black',label='Theory')

    try:

        Tpara = np.load(path_data+'A_0.3_Tpara_mu_bulk_D1.0_CFL'+CFL[i]+'.npy')
        Tperp = np.load(path_data+'A_0.3_Tperp_mu_bulk_D1.0_CFL'+CFL[i]+'.npy')
           
        Asim = Tpara/Tperp
        tsim = np.arange(tmin,len(Asim)*dt,dt)

        axs[i].plot(tsim,Asim,color='red',linestyle = '--',label='Simulation')
   
    except IOError:
        print('No file found')
        continue

    if i == 0:
        axs[i].legend(loc='best',bbox_to_anchor=(0.0,1.1,0.5,0.5),ncol=2)

    axs[i].text(8.8,0.75,panels[i]+' '+var+' = '+CFL[i],fontsize=12)

plt.savefig(path_fig+'A'+"{:.2f}".format(Ainit)+'vsSim.png',dpi=300)
print(path_fig+'A'+"{:.2f}".format(Ainit)+'vsSim.png')

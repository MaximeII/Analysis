import pytools as pt
import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

numproc = 20

runs = ['300_T2_B1','300_T3_B1','300_T4_B1','300_T5_B1','300_T2_B2','300_T3_B2','300_T4_B2','300_T5_B2','300_T2_B3','300_T3_B3','300_T4_B3','300_T5_B3']

loc = ['mid','right','left','top_mid','bot_mid','top_r','bot_r','top_l','bot_l']

path_fig  = '/wrk/users/dubart/analysis/subgrid/300_test/fig/'

fig, axs = plt.subplots(2,1,sharex=True)

legend_elements = [Line2D([0],[0],color='#1f77b4',lw=3,label='$T_{\\bot}/T_{\\parallel} = 2$'),
                   Line2D([0],[0],color='#ff7f0e',lw=3,label='$T_{\\bot}/T_{\\parallel} = 3$'),
                   Line2D([0],[0],color='#2ca02c',lw=3,label='$T_{\\bot}/T_{\\parallel} = 4$'),
                   Line2D([0],[0],color='#d62728',lw=3,label='$T_{\\bot}/T_{\\parallel} = 5$'),
                   Line2D([0],[0],color='black',lw=3,linestyle='solid',label='$\\beta_{\\parallel} = 1$'),
                   Line2D([0],[0],color='black',lw=3,linestyle='dotted',label='$\\beta_{\\parallel} = 2$'),
                   Line2D([0],[0],color='black',lw=3,linestyle=(0,(5,1)),label='$\\beta_{\\parallel} = 3$')]

axs[0].set_ylabel('$\\beta_{\\parallel}$',fontsize = 15)
axs[0].legend(handles=legend_elements,bbox_to_anchor=(-0.1,1.0,1.2,0.2),loc='upper right',ncol = 7,mode = 'expand',fontsize = 7)
axs[0].set_ylim(1.0,8.0)

axs[1].set_ylabel('$T_{\\bot}/T_{\\parallel}$',fontsize = 15)
axs[1].set_xlabel('Time (s)')
axs[1].set_xlim(0.0,300.0)
axs[1].set_ylim(1.0,5.0)

for i in range(0,len(runs)): 

    if i <= 3:
        linestyle = 'solid'
    elif i > 3 and i <= 7:
        linestyle = 'dotted'
    else:
        linestyle = (0,(5,1)) 

    if i == 0 or i == 4 or i == 8:
        color = '#1f77b4'
    elif i == 1 or i == 5 or i == 9:
        color = '#ff7f0e'
    elif i == 2 or i == 6 or i == 10:
        color = '#2ca02c'
    else:
        color = '#d62728'

    try:

        path_save = '/wrk/users/dubart/analysis/subgrid/300_test/data/'+runs[i]+'/'

        beta_par = np.load(path_save+'beta_par.npy')
        Taniso   = np.load(path_save+'T_aniso.npy')

        beta_avg   = np.mean(beta_par,axis=1)
        Taniso_avg = np.mean(Taniso,axis=1)
        
        time = np.linspace(0.0,len(beta_avg)/2,len(beta_avg))

        axs[0].plot(time,beta_avg,linestyle=linestyle,lw=2,color=color)
        axs[1].plot(time,Taniso_avg,linestyle=linestyle,lw=2,color=color)

    except:

        continue


plt.savefig(path_fig+'Variables_avg.png',dpi=300)
print(path_fig+'Variables_avg.png')


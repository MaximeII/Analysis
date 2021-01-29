import pytools as pt
import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

runs = ['300_T2_B1','300_T3_B1','300_T4_B1','300_T5_B1','300_T2_B2','300_T3_B2','300_T4_B2','300_T5_B2','300_T2_B3','300_T3_B3','300_T4_B3','300_T5_B3']

loc = ['mid','right','left','top_mid','bot_mid','top_r','bot_r','top_l','bot_l']

path_fig  = '/wrk/users/dubart/analysis/subgrid/300_test/fig/'

fig, ax = plt.subplots()

legend_elements = [Line2D([0],[0],color='#1f77b4',lw=3,label='$T_{\\bot}/T_{\\parallel} = 2$'),
                   Line2D([0],[0],color='#ff7f0e',lw=3,label='$T_{\\bot}/T_{\\parallel} = 3$'),
                   Line2D([0],[0],color='#2ca02c',lw=3,label='$T_{\\bot}/T_{\\parallel} = 4$'),
                   Line2D([0],[0],color='#d62728',lw=3,label='$T_{\\bot}/T_{\\parallel} = 5$'),
                   Line2D([0],[0],color='black',lw=3,linestyle='solid',label='$\\beta_{\\parallel} = 1$'),
                   Line2D([0],[0],color='black',lw=3,linestyle='dotted',label='$\\beta_{\\parallel} = 2$'),
                   Line2D([0],[0],color='black',lw=3,linestyle=(0,(5,1)),label='$\\beta_{\\parallel} = 3$')]

plt.ylabel('$\\sigma-\\sigma_0$',fontsize = 15)
plt.legend(handles=legend_elements,bbox_to_anchor=(-0.1,0.95,1.2,0.2),loc='upper right',ncol = 7,mode = 'expand',fontsize = 7)
plt.xlabel('Time (s)')
plt.xlim(0.5,50.0)

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


    path_save = '/wrk/users/dubart/analysis/subgrid/300_test/data/'+runs[i]+'/'

    std = np.load(path_save+'std_'+runs[i]+'.npy')
 
    time = np.linspace(0.5,50,len(std))

    plt.plot(time,std-std[0],linestyle=linestyle,lw=2,color=color)

    ax.xaxis.set_minor_locator(AutoMinorLocator())

plt.text(9.0,0.07,'$d\\sigma/dt \\approx 0.016$')
plt.text(18.0,0.10,'$d\\sigma/dt \\approx 0.010$')
plt.text(18.0,0.08,'$d\\sigma/dt \\approx 0.0085$')
plt.text(24.0,0.04,'$d\\sigma/dt \\approx 0.0055$')
plt.text(25.0,0.015,'$d\\sigma/dt \\approx 0.002$')

plt.savefig(path_fig+'std.png',dpi=300)
print(path_fig+'std.png')



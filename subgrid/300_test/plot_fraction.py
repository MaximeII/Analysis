import pytools as pt
import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

runs = ['300_T2_B1','300_T3_B1','300_T4_B1','300_T5_B1','300_T2_B2','300_T3_B2','300_T4_B2','300_T5_B2','300_T2_B3','300_T3_B3','300_T4_B3','300_T5_B3']

loc = ['mid','right','left','top_mid','bot_mid','top_r','bot_r','top_l','bot_l']

path_fig  = '/wrk/users/dubart/analysis/subgrid/300_test/fig/'

fig = plt.figure()

legend_elements = [Line2D([0],[0],color='#1f77b4',lw=3,label='$T_{\\bot}/T_{\\parallel} = 2$'),
                   Line2D([0],[0],color='#ff7f0e',lw=3,label='$T_{\\bot}/T_{\\parallel} = 3$'),
                   Line2D([0],[0],color='#2ca02c',lw=3,label='$T_{\\bot}/T_{\\parallel} = 4$'),
                   Line2D([0],[0],color='#d62728',lw=3,label='$T_{\\bot}/T_{\\parallel} = 5$'),
                   Line2D([0],[0],color='black',lw=3,linestyle='solid',label='$\\beta_{\\parallel} = 1$'),
                   Line2D([0],[0],color='black',lw=3,linestyle='dotted',label='$\\beta_{\\parallel} = 2$'),
                   Line2D([0],[0],color='black',lw=3,linestyle=(0,(5,1)),label='$\\beta_{\\parallel} = 1$')]

plt.ylabel('$\\langle \\mathrm{Fraction~of~phase~space~density} \\rangle$',fontsize = 15)
plt.legend(handles=legend_elements,bbox_to_anchor=(-0.1,0.95,1.2,0.2),loc='upper right',ncol = 7,mode = 'expand',fontsize = 7)
plt.xlabel('$\\mu = \\cos\\theta$')
plt.xlim(-1.0,1.0)
plt.ticklabel_format(axis='y',style='sci',scilimits=(-2,-2))


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
 
    lines = []
    
    filePA = open(path_save+'PA/mu'+str(1).rjust(7,'0')+'_0','r')
    
    for line in filePA:
        lines.append(line)
    
    angles   = np.fromstring(lines[3],dtype=float,sep=' ')
    binvalue = np.fromstring(lines[5],dtype=float,sep=' ')
    
    filePA.close()
    
    angles_center = (angles[0:-1]+angles[1:])/2

    Diff_fmu = np.load(path_save+'Diff_fmu_'+runs[i]+'.npy')
    
    plt.plot(angles_center,Diff_fmu,linestyle=linestyle,lw=2,color=color)


plt.savefig(path_fig+'Diff_fmu.png',dpi=300)
print(path_fig+'Diff_fmu.png')

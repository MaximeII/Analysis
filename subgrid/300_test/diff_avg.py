import pytools as pt
import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from multiprocessing import Pool

# Avoids opening a figure window
if str(matplotlib.get_backend()) != 'Agg':
    plt.switch_backend('Agg')

numproc = 20

runs = ['300_T2_B1','300_T3_B1','300_T4_B1','300_T5_B1','300_T2_B2','300_T3_B2','300_T4_B2','300_T5_B2','300_T2_B3','300_T3_B3','300_T4_B3','300_T5_B3']

loc = ['mid','right','left','top_mid','bot_mid','top_r','bot_r','top_l','bot_l']

path_fig  = '/wrk/users/dubart/analysis/subgrid/300_test/fig/'

bulkStart = 1
bulkEnd   = 600

def plot_avg_diff(step):

    fig, axs = plt.subplots(2,1,sharex=True)
    
    legend_elements = [Line2D([0],[0],color='#1f77b4',lw=3,label='$T_{\\bot}/T_{\\parallel} = 2$'),
                       Line2D([0],[0],color='#ff7f0e',lw=3,label='$T_{\\bot}/T_{\\parallel} = 3$'),
                       Line2D([0],[0],color='#2ca02c',lw=3,label='$T_{\\bot}/T_{\\parallel} = 4$'),
                       Line2D([0],[0],color='#d62728',lw=3,label='$T_{\\bot}/T_{\\parallel} = 5$'),
                       Line2D([0],[0],color='black',lw=3,linestyle='solid',label='$\\beta_{\\parallel} = 1$'),
                       Line2D([0],[0],color='black',lw=3,linestyle='dotted',label='$\\beta_{\\parallel} = 2$'),
                       Line2D([0],[0],color='black',lw=3,linestyle=(0,(3,1,1,1,1,1)),label='$\\beta_{\\parallel} = 1$')]
    
    axs[0].set_ylabel('$f(\\mu)$',fontsize = 15)
    axs[0].legend(handles=legend_elements,bbox_to_anchor=(-0.1,1.0,1.2,0.2),loc='upper right',ncol = 7,mode = 'expand',fontsize = 7)
    axs[0].set_ylim(1.0e3,3.0e5)
    axs[0].ticklabel_format(axis='y',style='sci',scilimits=(4,5))
    
    axs[1].set_ylabel('$D_{\\mu\\mu}$',fontsize = 15)
    axs[1].set_xlabel('$\\mu = \\cos{\\theta}$',fontsize=15)
    axs[1].set_xlim(-1.0,1.0)
    axs[1].set_ylim(0.0,0.10)
    
    fig.suptitle('t = '+str(step/2.0)+' s',fontsize = 15)
    
    for i in range(0,len(runs)): 
    
        if i <= 3:
            linestyle = 'solid'
        elif i > 3 and i <= 7:
            linestyle = 'dotted'
        else:
            linestyle = (0,(3,1,1,1,1,1)) 
    
        try:
    
            path_save = '/wrk/users/dubart/analysis/subgrid/300_test/data/'+runs[i]+'/'
     
            # f(mu)
            lines = []

            filePA = open(path_save+'PA/mu'+str(step).rjust(7,'0')+'_'+str(i),'r')

            for line in filePA:
                lines.append(line)

            angles   = np.fromstring(lines[3],dtype=float,sep=' ')
            binvalue = np.fromstring(lines[5],dtype=float,sep=' ')

            filePA.close()

            angles_center = (angles[0:-1]+angles[1:])/2

            
            binvalue_tmp = np.zeros(len(binvalue))
            
            for j in range(0,9):
                lines = []

                filePA = open(path_save+'PA/mu'+str(step).rjust(7,'0')+'_'+str(j),'r')

                for line in filePA:
                    lines.append(line)

                binvalue = np.fromstring(lines[5],dtype=float,sep=' ')

                filePA.close()

                binvalue_tmp = binvalue_tmp + binvalue

            binvalue_avg = binvalue_tmp / 9.0
  
            np.save(path_save+'fmu_avg.npy',binvalue_avg)
      
            axs[0].plot(angles_center,binvalue_avg,linestyle=linestyle,lw=2)

            Dmumu_tmp = np.zeros(len(binvalue))            

            for j in range(0,9):

                Dmumu = np.load(path_save+'Dmumu_'+str(j)+'.npy')

                Dmumu = Dmumu[step-bulkStart,:]

                Dmumu_tmp = Dmumu_tmp + Dmumu

            Dmumu_avg = Dmumu_tmp / 9.0
    
            axs[1].plot(angles_center,Dmumu_avg,linestyle=linestyle,lw=2)
    
        except:
    
            continue
    
    
    plt.savefig(path_fig+'Dmumu_avg_'+str(step).rjust(7,'0')+'.png',dpi=300)
    print(path_fig+'Dmumu_avg_'+str(step).rjust(7,'0')+'.png')
    plt.close()

pool = Pool(numproc)
pool.map(plot_avg_diff,range(bulkStart,bulkEnd+1))
pool.terminate()

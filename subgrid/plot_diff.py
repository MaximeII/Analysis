import pytools as pt
import numpy as np
import sys
from multiprocessing import Pool
import matplotlib
import matplotlib.pyplot as plt

# Avoids opening a figure window
if str(matplotlib.get_backend()) != 'Agg':
    plt.switch_backend('Agg')


numproc = 20 

run = sys.argv[1]
bulkStart = int(sys.argv[2])
bulkEnd   = int(sys.argv[3])

#path_bulk = '/wrk/users/dubart/'+run+'km/Shock/'
#path_bulk = '/wrk/users/dubart/'+run+'km/maxwellian_periodic/'
path_bulk = '/wrk/users/dubart/'+run+'km/bulk/'
path_save = '/wrk/users/dubart/analysis/subgrid/data/'+run+'/periodic/'
path_fig  = '/wrk/users/dubart/analysis/subgrid/fig/'+run+'/periodic/'

Dmumu = np.load(path_save+'Dmumu.npy')

def plot_diff(step):

    lines = []

    filePA = open(path_save+'PA/mu'+str(step).rjust(7,'0'),'r')

    for line in filePA:
        lines.append(line)

    angles   = np.fromstring(lines[3],dtype=float,sep=' ')
    binvalue = np.fromstring(lines[5],dtype=float,sep=' ')

    filePA.close()

    angles_center = (angles[0:-1]+angles[1:])/2

    fig, axs = plt.subplots(2,1,sharex=True)

    axs[0].plot(angles_center,binvalue,linewidth=2)
    axs[0].set_ylabel('Density',fontsize = 15)
    axs[0].set_ylim(1e3,1e5)
    axs[0].ticklabel_format(axis='y',style='sci',scilimits=(4,5))


    axs[1].plot(angles_center,abs(Dmumu[step - bulkStart,:]),linewidth=2)
    axs[1].set_xlabel('$\\mu = \\cos{\\theta}$',fontsize = 15)
    axs[1].set_ylabel('$D_{\\mu \\mu}$',fontsize = 15)
    axs[1].set_xlim(-1.0,1.0)
    axs[1].set_ylim(-0.10,0.10)

    fig.suptitle('t = '+str(step/2.0)+' s',fontsize = 15)

    plt.savefig(path_fig+'PA_DC_'+str(step).rjust(7,'0')+'.png',dpi=300)
    print(path_fig+'PA_DC_'+str(step).rjust(7,'0')+'.png')
    plt.close()

    return

pool = Pool(numproc)
pool.map(plot_diff,range(bulkStart,bulkEnd+1))
pool.terminate()
print('Job Done')    

import pytools as pt
import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

run = sys.argv[1]

bulkStart = int(sys.argv[2])
bulkEnd   = int(sys.argv[3])

path_save = '/wrk/users/dubart/analysis/subgrid/300_test/Green_diff/data/'
path_fig  = '/wrk/users/dubart/analysis/subgrid/300_test/Green_diff/fig/diffusion/restart/'
path_bulk = '/wrk/users/dubart/300_test/900km/diffusion/restart/bulk/'

RE = 6371e3 # m

fig, axs = plt.subplots(3,1,sharex=True)

rho = np.load(path_save+'rho_'+run+'.npy')

beta      = np.load(path_save+'beta_'+run+'.npy')
beta_para = np.load(path_save+'beta_para_'+run+'.npy') 
beta_perp = np.load(path_save+'beta_perp_'+run+'.npy')

Temperature = np.load(path_save+'T_'+run+'.npy')
T_para      = np.load(path_save+'Tpara_'+run+'.npy')    
T_perp      = np.load(path_save+'Tperp_'+run+'.npy')
T_aniso     = T_para/T_perp

time = np.linspace(bulkStart/2,bulkEnd/2,len(beta_para))

axs[0].plot(time,rho,linewidth = 2)
axs[0].set_ylabel('n [m$^{-3}$]',fontsize=15)
axs[0].ticklabel_format(axis='y',style="sci",scilimits=(6,6))
axs[0].yaxis.set_minor_locator(AutoMinorLocator())

axs[1].plot(time,T_aniso,linewidth = 2)
axs[1].set_ylabel('$T_\\parallel/T_\\bot$',fontsize=15)
axs[1].yaxis.set_minor_locator(AutoMinorLocator())

axs[2].plot(time,beta,linewidth=2)
axs[2].set_ylabel('$\\beta$',fontsize=15)
axs[2].yaxis.set_minor_locator(AutoMinorLocator())

axs[2].set_xlabel('Time (s)')
axs[2].set_xlim(bulkStart/2,bulkEnd/2)

plt.savefig(path_fig+'Variables_'+run+'.png',dpi=300)
print(path_fig+'Variables_'+run+'.png')

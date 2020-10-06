import numpy as np
import matplotlib.pyplot as plt
import sys

theta = 89.2

path_output = '/wrk/users/dubart/analysis/hydros/Lucile/output/'
path_fig    = '/wrk/users/dubart/analysis/hydros/Lucile/fig/'

HYDROS=np.loadtxt(path_output+'Lucile_1D_'+str(theta)+'.log',dtype=float)

k     = HYDROS[:,0]
omega = HYDROS[:,2]

plt.plot(k,omega,linewidth=2)
plt.xlabel('$k*d_i$',fontsize=15)
plt.ylabel('$\\omega/\\Omega_c$',fontsize=15)
plt.ylim(0.0,4.0)
plt.xlim(0.0,5.0)

plt.savefig(path_fig+'89.2_0.01.png',dpi=300)
print(path_fig+'89.2_0.01.png')

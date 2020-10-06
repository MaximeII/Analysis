import numpy as np
import matplotlib.pyplot as plt
import sys
from mpl_toolkits import mplot3d

theta = np.arange(40.0,90.1,0.1)
k     = np.arange(0.01,5.1,0.01)

w_start = 0.01

path_output = '/wrk/users/dubart/analysis/hydros/Lucile/output/'
path_fig    = '/wrk/users/dubart/analysis/hydros/Lucile/fig/'

fig = plt.figure()
ax  = plt.axes(projection='3d')
ax.set_xlabel('$k * d_i$')
ax.set_ylabel('$\\theta$')
ax.set_zlabel('$\\omega / \\Omega_C$')

for i in range(0,len(theta)):

    HYDROS=np.loadtxt(path_output+'Lucile_'+str(w_start)+'_'+str(round(theta[i],1))+'.log',dtype=float)
    print('Lucile_'+str(w_start)+'_'+str(round(theta[i],1))+'.log') 
   
    omega    = np.empty(len(k))
    omega[:] = np.NAN   

    if len(HYDROS.shape) == 2:
        omegaHY = HYDROS[:,2]
        omega[:len(omegaHY)] = omegaHY
    else:
        omegaHY = HYDROS[2]
        omega[0] = omegaHY

    ax.scatter3D(k,theta[i],omega)

plt.savefig(path_fig+'Lucile_3D_'+str(w_start)+'.png',DPI=300)  
print(path_fig+'Lucile_3D_'+str(w_start)+'.png')
    

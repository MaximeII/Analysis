import numpy as np
import matplotlib.pyplot as plt
import sys
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D

theta     = np.arange(40.0,90.1,0.1)
k         = np.arange(0.00,5.01,0.01)
omegaFULL = np.empty([len(k),len(theta)])
gammaFULL = np.empty([len(k),len(theta)])

w_start = 1.0

path_Owen = '/wrk/users/dubart/analysis/hydros/Lucile/Owen/'
path_output = '/wrk/users/dubart/analysis/hydros/Lucile/output/'
path_fig    = '/wrk/users/dubart/analysis/hydros/Lucile/fig/'

fig = plt.figure()
ax  = fig.add_subplot(111,projection='3d')
ax.set_xlabel('$k * d_i$')
ax.set_ylabel('$\\theta$')
ax.set_zlabel('$\\omega / \\Omega_C$')
ax.set_zlim(0.0,2.0)

for i in range(0,len(theta)):

    HYDROS=np.loadtxt(path_output+'Lucile_'+str(w_start)+'_'+str(round(theta[i],1))+'.log',dtype=float)
#    print('Lucile_'+str(w_start)+'_'+str(round(theta[i],1))+'.log') 
   
    omega    = np.empty(len(k))
    omega[:] = np.NAN   
    
    gamma    = np.empty(len(k))
    gamma[:] = np.NAN

    if len(HYDROS.shape) == 2:
        omegaHY = HYDROS[:,2]
        omega[:len(omegaHY)] = omegaHY
        gammaHY = -HYDROS[:,3]
        gamma[:len(gammaHY)] = gammaHY
    else:
        omegaHY  = HYDROS[2]
        omega[0] = omegaHY
        gammaHY  = -HYDROS[3]
        gamma[0] = gammaHY

    gammaFULL[:,i] = gamma
    omegaFULL[:,i] = omega

masked_omega=np.ma.masked_invalid(omegaFULL)
masked_gamma=np.ma.masked_invalid(gammaFULL)

w_min = np.amin(masked_omega)
w_max = 2.0

gamma_min = 0.0
gamma_max = 0.5

for i in range(0,len(theta)):

    HYDROS=np.loadtxt(path_output+'Lucile_'+str(w_start)+'_'+str(round(theta[i],1))+'.log',dtype=float)
#    print('Lucile_'+str(w_start)+'_'+str(round(theta[i],1))+'.log') 

    omega    = np.empty(len(k))
    omega[:] = np.NAN

    if len(HYDROS.shape) == 2:
        omegaHY = HYDROS[:,2]
        omega[:len(omegaHY)] = omegaHY
        gammaHY = -HYDROS[:,3]
        gamma[:len(gammaHY)] = gammaHY
    else:
        omegaHY  = HYDROS[2]
        omega[0] = omegaHY
        gammaHY  = -HYDROS[3]
        gamma[0] = gammaHY
    
    p=ax.scatter3D(k,theta[i],omega,c=gamma,cmap='viridis',s=2,vmin=gamma_min,vmax=gamma_max)

cbar=fig.colorbar(p)
#cbar.mappable.set_clim(0.0,2.0)
cbar.set_label('$\\gamma/\\Omega_C$')

Owen1  = np.loadtxt(path_Owen+'2019 01 11 03 25 002019 01 11 03 31 59.txt',skiprows=1,dtype=float)
omega  = Owen1[:,1]
kowen  = Owen1[:,2]
thetaO = Owen1[:,3]

for i in range(0,len(thetaO)):
    if thetaO[i] > 90.0:
        thetaO[i] = abs(thetaO[i] - 180.0)

ax.scatter3D(kowen,thetaO,omega,s=1,color='red',marker='x')
ax.set_ylim(40.0,90.0)

plt.savefig(path_fig+'Lucile_3D_'+str(w_start)+'.png',DPI=300)  
print(path_fig+'Lucile_3D_'+str(w_start)+'.png')
    

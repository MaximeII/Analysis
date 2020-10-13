import numpy as np
import matplotlib.pyplot as plt
import sys
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D

path_Owen = '/wrk/users/dubart/analysis/hydros/Lucile/Owen/'
path_output = '/wrk/users/dubart/analysis/hydros/Lucile/output/'
path_fig    = '/wrk/users/dubart/analysis/hydros/Lucile/fig/'

Owen1  = np.loadtxt(path_Owen+'2019 01 11 03 25 002019 01 11 03 31 59.txt',skiprows=1,dtype=float)
omega  = Owen1[:,1]
kowen  = Owen1[:,2]
thetaO = Owen1[:,3]

Owen2  = np.loadtxt(path_Owen+'2019 01 11 03 40 002019 01 11 03 44 59.txt',skiprows=1,dtype=float)
omega  = np.append(omega,Owen2[:,1])
kowen  = np.append(kowen,Owen2[:,2])
thetaO = np.append(thetaO,Owen2[:,3])

Owen3  = np.loadtxt(path_Owen+'2019 01 11 03 50 002019 01 11 03 51 59.txt',skiprows=1,dtype=float)
omega  = np.append(omega,Owen3[:,1])
kowen  = np.append(kowen,Owen3[:,2])
thetaO = np.append(thetaO,Owen3[:,3])

omega = abs(omega)

for i in range(0,len(thetaO)):
    if thetaO[i] > 90.0:
        thetaO[i] = abs(thetaO[i] - 180.0)

theta = np.arange(5.0,90.0,5.0)

for i in range(0,len(theta)):
    
    fig,axs = plt.subplot(2,1)
    
    w_start = 1.0

    HYDROS=np.loadtxt(path_output+'Lucile_'+str(w_start)+'_'+str(round(theta[i],1))+'.log',dtype=float)

    kHY     = HYDROS[:,1]
    omegaHY = HYDROS[:,2]
    gammaHY = HYDROS[:,3]

    axs[0].plot(kHY,omegaHY,color='blue',linewidth=2)
    axs[0].plot(kHY,gammaHY,color='blue',linewidth=2,linestyle='--')
    axs[0].set_xlim(0.0,3.0)
    axs[0].set_ylim(0.0,2.0)
    axs[0].set_xlabel('$k*d_i$')
    axs[0].set_ylabel('$\\omega/\\Omega_C$')
    axs[0].set_title('\\omega_{\\mathrm{guess}} = 1.0')

    w_start = 0.01

    HYDROS=np.loadtxt(path_output+'Lucile_'+str(w_start)+'_'+str(round(theta[i],1))+'.log',dtype=float)

    kHY     = HYDROS[:,1]
    omegaHY = HYDROS[:,2]
    gammaHY = HYDROS[:,3]

    axs[1].plot(kHY,omegaHY,color='orange',linewidth=2)
    axs[1].plot(kHY,gammaHY,color='orange',linewidth=2,linestyle='--')
    axs[1].set_xlim(0.0,3.0)
    axs[1].set_ylim(0.0,2.0)
    axs[1].set_xlabel('$k*d_i$')
    axs[1].set_ylabel('$\\omega/\\Omega_C$')
    axs[1].set_title('\\omega_{\\mathrm{guess}} = 0.01')
    
    #thetapts = thetaO[(thetaO>=(theta[i]-5.0)) & (thetaO<=(theta[i]+5.0))]
    kpts     = kowen[(thetaO>=(theta[i]-5.0)) & (thetaO<=(theta[i]+5.0))]
    omegapts = omega[(thetaO>=(theta[i]-5.0)) & (thetaO<=(theta[i]+5.0))]

    axs[0].scatter(kpts,omegapts,color='black',marker='x')
    axs[1].scatter(kpts,omegapts,color='black',marker='x')

    plt.savefig(path_fig+'Lucile_2D_'+str(round(theta[i],1))+'.png',DPI=300)
    print(path_fig+'Lucile_2D_'+str(round(theta[i],1))+'.png')

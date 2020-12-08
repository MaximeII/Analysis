import numpy as np
import matplotlib.pyplot as plt
import sys
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D

freq = sys.argv[1]

path_Owen = '/wrk/users/dubart/analysis/hydros/Lucile/Owen/'
path_output = '/wrk/users/dubart/analysis/hydros/Lucile/output/'
path_fig    = '/wrk/users/dubart/analysis/hydros/Lucile/fig/'

Owen1   = np.loadtxt(path_Owen+'2019 01 11 03 25 002019 01 11 03 31 59.txt',skiprows=1,dtype=float)
f_sc1   = Owen1[:,0]
omega1  = abs(Owen1[:,1])
kowen1  = Owen1[:,2]
thetaO1 = Owen1[:,3]

Owen2   = np.loadtxt(path_Owen+'2019 01 11 03 40 002019 01 11 03 44 59.txt',skiprows=1,dtype=float)
f_sc2   = Owen2[:,0]
omega2  = abs(Owen2[:,1])
kowen2  = Owen2[:,2]
thetaO2 = Owen2[:,3]

Owen3   = np.loadtxt(path_Owen+'2019 01 11 03 50 002019 01 11 03 51 59.txt',skiprows=1,dtype=float)
f_sc3   = Owen3[:,0]
omega3  = abs(Owen3[:,1])
kowen3  = Owen3[:,2]
thetaO3 = Owen3[:,3]

for i in range(0,len(thetaO1)):
    if thetaO1[i] > 90.0:
        thetaO1[i] = abs(thetaO1[i] - 180.0)

for i in range(0,len(thetaO2)):
    if thetaO2[i] > 90.0:
        thetaO2[i] = abs(thetaO2[i] - 180.0)

for i in range(0,len(thetaO3)):
    if thetaO3[i] > 90.0:
        thetaO3[i] = abs(thetaO3[i] - 180.0)

if freq == 'low':

    omega1  = omega1[f_sc1 <= 0.1]  
    kowen1  = kowen1[f_sc1 <= 0.1]
    thetaO1 = thetaO1[f_sc1 <= 0.1]

    omega2  = omega2[f_sc2 <= 0.1]
    kowen2  = kowen2[f_sc2 <= 0.1]
    thetaO2 = thetaO2[f_sc2 <= 0.1]

    omega3  = omega3[f_sc3 <= 0.1]
    kowen3  = kowen3[f_sc3 <= 0.1]
    thetaO3 = thetaO3[f_sc3 <= 0.1]

    sign = '<='

    xmax = 0.5

elif freq == 'high':

    omega1  = omega1[f_sc1 > 0.1]
    kowen1  = kowen1[f_sc1 > 0.1]
    thetaO1 = thetaO1[f_sc1 > 0.1]

    omega2  = omega2[f_sc2 > 0.1]
    kowen2  = kowen2[f_sc2 > 0.1]
    thetaO2 = thetaO2[f_sc2 > 0.1]

    omega3  = omega3[f_sc3 > 0.1]
    kowen3  = kowen3[f_sc3 > 0.1]
    thetaO3 = thetaO3[f_sc3 > 0.1]    

    sign = '>'

    xmax = 3.0

theta = np.arange(1.5,90.0,3.0)

for i in range(0,len(theta)):
    
    fig,axs = plt.subplots(2,1)
    plt.subplots_adjust(hspace=0.5)

    w_start = 1.0

    HYDROS=np.loadtxt(path_output+'Lucile_'+str(w_start)+'_'+str(round(theta[i],1))+'.log',dtype=float)

    kHY     = HYDROS[:,0]
    omegaHY = HYDROS[:,2]
    gammaHY = - HYDROS[:,3]

    axs[0].plot(kHY,omegaHY,color='blue',linewidth=2,label='$\\omega$')
    axs[0].plot(kHY,gammaHY,color='blue',linewidth=2,linestyle='--',label='$\\gamma$')
    axs[0].set_xlim(0.0,xmax)
    axs[0].set_ylim(0.0,2.0)
    axs[0].set_xlabel('$k*d_i$')
    axs[0].set_ylabel('$\\omega/\\Omega_C$')
    axs[0].set_title('$\\omega_{\\mathrm{guess}} = 1.0$, $\\theta = '+str(round(theta[i],1))+'$')

    w_start = 0.01

    HYDROS=np.loadtxt(path_output+'Lucile_'+str(w_start)+'_'+str(round(theta[i],1))+'.log',dtype=float)

    if len(HYDROS.shape) == 2:
        kHY     = HYDROS[:,0]
        omegaHY = HYDROS[:,2]
        gammaHY = - HYDROS[:,3]
    else:
        kHY     = HYDROS[0]
        omegaHY = HYDROS[2]
        gammaHY = - HYDROS[3]

    axs[1].plot(kHY,omegaHY,color='orange',linewidth=2)
    axs[1].plot(kHY,gammaHY,color='orange',linewidth=2,linestyle='--')
    axs[1].set_xlim(0.0,xmax)
    axs[1].set_ylim(0.0,2.0)
    axs[1].set_xlabel('$k*d_i$')
    axs[1].set_ylabel('$\\omega/\\Omega_C$')
    axs[1].set_title('$\\omega_{\\mathrm{guess}} = 0.01$')
    
    kpts1     = kowen1[(thetaO1>=(theta[i]-1.5)) & (thetaO1<=(theta[i]+1.5))]
    omegapts1 = omega1[(thetaO1>=(theta[i]-1.5)) & (thetaO1<=(theta[i]+1.5))]

    kpts2     = kowen2[(thetaO2>=(theta[i]-1.5)) & (thetaO2<=(theta[i]+1.5))]
    omegapts2 = omega2[(thetaO2>=(theta[i]-1.5)) & (thetaO2<=(theta[i]+1.5))]

    kpts3     = kowen3[(thetaO3>=(theta[i]-1.5)) & (thetaO3<=(theta[i]+1.5))]
    omegapts3 = omega3[(thetaO3>=(theta[i]-1.5)) & (thetaO3<=(theta[i]+1.5))]
    
    axs[0].scatter(kpts1,omegapts1,color='black',marker='x',label='25 to 31')
    axs[1].scatter(kpts1,omegapts1,color='black',marker='x')

    axs[0].scatter(kpts2,omegapts2,color='black',marker='o',label='40 to 44')
    axs[1].scatter(kpts2,omegapts2,color='black',marker='o')
    
    axs[0].scatter(kpts3,omegapts3,color='black',marker='^',label='50 to 51')
    axs[1].scatter(kpts3,omegapts3,color='black',marker='^')
    
    axs[0].legend(loc='best')

    plt.suptitle("$f_{sc} "+sign+" 0.1~\\mathrm{Hz}$")

    plt.savefig(path_fig+'Lucile_2D_'+str(round(theta[i],1))+'_'+freq+'.png',DPI=300)
    print(path_fig+'Lucile_2D_'+str(round(theta[i],1))+'_'+freq+'.png')

    plt.close()

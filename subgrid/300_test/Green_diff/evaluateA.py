import numpy as np
import scipy.integrate
import matplotlib
import matplotlib.pyplot as plt
import sys
import pytools as pt

path_fig  = '/wrk/users/dubart/analysis/subgrid/300_test/Green_diff/fig/evaluateA/'
path_save = '/wrk/users/dubart/analysis/subgrid/300_test/Green_diff/data/'
path_bulk = '/wrk/users/dubart/diff_test/proc_test/mu_bulk/'

murange  = np.linspace(-1.0,1.0,101)
dmu      = abs(murange[0] - murange[1])

Dmumu = float(sys.argv[1])

f = pt.vlsvfile.VlsvReader(path_bulk+'initial-grid.0000000.vlsv')

Tpara_init = f.read_variable('vg_t_parallel')
Tperp_init = f.read_variable('vg_t_perpendicular')

n_init     = 3.0e6
m          = 1.6726219e-27
kB         = 1.38064852e-23

dt   = float(sys.argv[2])
time = np.arange(0.0,12.0,dt)

A    = np.zeros(len(time))
A[0] = Tpara_init/Tperp_init

vrange = np.linspace(-10.0,10.0,201)*1e6
dv     = abs(vrange[0] - vrange[1])

v2D      = np.sqrt(2.0*vrange**2)
vmurange = np.linspace(min(v2D),max(v2D),len(vrange))
dvmu     = abs(vmurange[0] - vmurange[1])
v,mu     = np.meshgrid(vmurange,murange)

fvmu = 1.0/np.sqrt(2.0*np.pi) * (m/kB)**(3.0/2.0) * (n_init*v**2)/(np.sqrt(Tpara_init)*Tperp_init) * np.exp(-m/(2.0*kB) * ((v**2 *mu**2)/Tpara_init + (v**2 * (1-mu**2))/Tperp_init))



for t in range(1,len(time)):

    dfdt    = np.zeros(fvmu.shape)
    dfdmu   = np.zeros(fvmu.shape)
    dfdmumu = np.zeros(fvmu.shape)
    
    dfdmu[1:-1,:]   = (fvmu[2:,:] - fvmu[:-2,:])/(2.0*dmu)
    dfdmu[0,:]      = (fvmu[1,:]  - fvmu[0,:])/dmu
    dfdmu[-1,:]     = (fvmu[-1,:] - fvmu[-2,:])/dmu
    
    for j in range(1,len(murange)-1):
        dfdmumu[j,:] = (fvmu[j+1,:] - 2.0*fvmu[j,:] + fvmu[j-1,:])/(dmu*dmu)

    dfdmumu[0,:]  = (fvmu[2,:] - 2.0*fvmu[1,:] + fvmu[0,:])/(dmu*dmu)
    dfdmumu[-1,:] = (fvmu[-3,:] - 2.0*fvmu[-2,:] + fvmu[-1,:])/(dmu*dmu)
    
    for j in range(0,len(murange)):
        dfdt[j,:] = Dmumu * ( -2*mu[j,:]*dfdmu[j,:] + (1.0 - mu[j,:]*mu[j,:])*dfdmumu[j,:])
        fvmu[j,:] = fvmu[j,:] + dfdt[j,:]*dt
    
    n = scipy.integrate.simps(scipy.integrate.simps(fvmu,dx=dvmu),dx=dmu)

    Tpara = (m/(n*kB)) * scipy.integrate.simps(scipy.integrate.simps(v**2 * mu**2 * fvmu,dx=dvmu),dx=dmu)
    Tperp = (m/(2.0*n*kB))* scipy.integrate.simps(scipy.integrate.simps(v**2 * (1-mu**2) * fvmu,dx=dvmu),dx=dmu)

    A[t] = Tpara/Tperp
     
np.save(path_save+'A_'+"{:.2f}".format(A[0])+'_D'+str(Dmumu)+'_dt'+str(dt)+'.npy',A)
print(path_save+'A_'+"{:.2f}".format(A[0])+'_D'+str(Dmumu)+'_dt'+str(dt)+'.npy')

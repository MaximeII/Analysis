import numpy as np
import scipy.integrate
import matplotlib
import matplotlib.pyplot as plt
import sys
import pytools as pt

run = sys.argv[1]

path_save = '/wrk/users/dubart/analysis/subgrid/300_test/Green_diff/data/mu/'+run+'/'
path_bulk = '/wrk/users/dubart/diff_test/dmumu/'+run+'/bulk/'

murange  = np.linspace(-1.0,1.0,101)
dmu      = abs(murange[0] - murange[1])

Dmumu = float(sys.argv[2])

CellID = int(np.load(path_save+'CellID_'+run+'.npy'))

f = pt.vlsvfile.VlsvReader(path_bulk+'initial-grid.0000000.vlsv')

Tpara_init = f.read_variable('vg_t_parallel',CellID)
Tperp_init = f.read_variable('vg_t_perpendicular',CellID)

n_init     = 3.0e6
m          = 1.6726219e-27
kB         = 1.38064852e-23

dt   = float(sys.argv[3])
time = np.arange(0.0,12.0,dt)

vrange = np.linspace(-10.0,10.0,201)*1e6
dv     = abs(vrange[0] - vrange[1])

v2D      = np.sqrt(2.0*vrange**2)
vmurange = np.linspace(min(v2D),max(v2D),len(vrange))
dvmu     = abs(vmurange[0] - vmurange[1])
v,mu     = np.meshgrid(vmurange,murange)

fvmu = 1.0/np.sqrt(2.0*np.pi) * (m/kB)**(3.0/2.0) * (n_init*v**2)/(np.sqrt(Tpara_init)*Tperp_init) * np.exp(-m/(2.0*kB) * ((v**2 *mu**2)/Tpara_init + (v**2 * (1-mu**2))/Tperp_init))

avg_fmu = np.mean(fvmu,axis=1)
avg_mu  = np.mean(mu,axis=1)

fmu = np.empty([len(avg_fmu),len(time)])
fmu[:,0] = avg_fmu
Dmumu_solve = np.empty([len(avg_fmu),len(time)])
Dmumu_solve[:,0] = (1-avg_mu**2)*Dmumu

for t in range(1,len(time)):

    dfdt    = np.zeros(avg_fmu.shape)
    dfdmu   = np.zeros(avg_fmu.shape)
    dfdmumu = np.zeros(avg_fmu.shape)
    
    dfdmu[1:-1]   = (avg_fmu[2:] - avg_fmu[:-2])/(2.0*dmu)
    dfdmu[0]      = (avg_fmu[1]  - avg_fmu[0])/dmu
    dfdmu[-1]     = (avg_fmu[-1] - avg_fmu[-2])/dmu
    
    # dfdmumu
    dfdmumu[0]    = (avg_fmu[2]  - 2*avg_fmu[1]    + avg_fmu[0]) / dmu**2
    dfdmumu[1:-1] = (avg_fmu[2:] - 2*avg_fmu[1:-1] + avg_fmu[:-2]) / dmu**2
    dfdmumu[-1]   = (avg_fmu[-1] - 2*avg_fmu[-2]   + avg_fmu[-3]) / dmu**2

    
    for j in range(0,len(murange)):
        dfdt[j]   = Dmumu * ( -2*avg_mu[j]*dfdmu[j] + (1.0 - avg_mu[j]**2)*dfdmumu[j])
        avg_fmu[j] = avg_fmu[j] + dfdt[j]*dt

    fmu[:,t] = avg_fmu

    Matrix = np.diag(dfdmumu,0) + np.diag(dfdmu[:-1]/(2*dmu),1) + np.diag(-dfdmu[1:]/(2*dmu),-1)
    try:
        InvM = np.linalg.inv(Matrix)
        Dmumu_solve[:,t] = np.matmul(InvM,dfdt)
    except Exception as e:
        print('t = '+str(t)+': '+str(e))

np.save(path_save+'fmu_1D_test.npy',fmu)
print(path_save+'fmu_1D_test.npy')
np.save(path_save+'Dmumu_1D_test.npy',Dmumu_solve)
print(path_save+'Dmumu_1D_test.npy')

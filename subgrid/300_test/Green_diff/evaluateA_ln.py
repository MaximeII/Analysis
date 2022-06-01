import numpy as np
import scipy.integrate
import matplotlib
import matplotlib.pyplot as plt
import sys

path_fig  = '/wrk/users/dubart/analysis/subgrid/300_test/Green_diff/fig/evaluateA/'
path_save = '/wrk/users/dubart/analysis/subgrid/300_test/Green_diff/data/'
 
v     = float(sys.argv[1])
dt    = float(sys.argv[2])
Dmumu = float(sys.argv[3])

mu  = np.linspace(-1.0,1.0,101)
dmu = abs(mu[0] - mu[1])
 
Ainit = float(sys.argv[4])

time  = np.arange(0.0,12.0,dt)
tplot = np.linspace(0.0,12.0,101)

A    = np.zeros(len(time))
A[0] = Ainit

#vrange = np.linspace(-10.0,10.0,len(mu))

#vpara,vperp = np.meshgrid(vrange,vrange)

f = 2*np.pi * v**2 *np.exp(- v*v * mu*mu - A[0] * (v*v * (1.0 - mu*mu)))

#fplot =  np.exp(- vpara*vpara - A[0] * vperp*vperp)

print(f)
print(Ainit)
#ax_fmu = plt.subplot(221)
#ax_fmu.plot(mu,f,color='black')
#ax_fmu.set_xlim(-1.0,1.0)
#ax_fmu.set_ylim(0.0,0.001)
#ax_fmu.set_xlabel('$\\mu$',fontsize=15)
#ax_fmu.set_ylabel('$f$',fontsize=15)
#
#ax_fv = plt.subplot(222)
#ax_fv.pcolormesh(vpara,vperp,fplot,cmap='viridis',shading='auto')
#ax_fv.set_xlim(vrange[0]/3,vrange[-1]/3)
#ax_fv.set_ylim(vrange[0]/3,vrange[-1]/3)
#ax_fv.set_xlabel('$v_\\parallel$',fontsize=15)
#ax_fv.set_ylabel('$v_\\bot$',fontsize=15)
#ax_fv.set_aspect('equal')
#
#ax_A = plt.subplot(212)
#ax_A.scatter(time[0],A[0],s=10)
#ax_A.set_xlim(0.0,12.0)
#ax_A.set_ylim(A[0],1.0)
#ax_A.set_xlabel('Time (s)',fontsize=15)
#ax_A.set_ylabel('A',fontsize=15)
#
#plt.subplots_adjust(hspace = 0.4)
#
#plt.savefig(path_fig+'A_v'+str(v)+'_D'+str(Dmumu)+'_t'+"{:.3f}".format(time[0])+'.png',dpi=300)
#print(path_fig+'A_v'+str(v)+'_D'+str(Dmumu)+'_t'+"{:.3f}".format(time[0])+'.png')
#plt.close()

for t in range(1,len(time)):

    dfdt    = np.zeros(len(mu))
    dfdmu   = np.zeros(len(mu))
    dfdmumu = np.zeros(len(mu))
    
    dfdmu[1:-1]   = (f[2:] - f[:-2])/(2*dmu)
    dfdmu[0]      = (f[1] - f[0])/dmu
    dfdmu[-1]     = (f[-1] - f[-2])/dmu
    
    for j in range(1,len(mu)-1):
        dfdmumu[j] = (f[j+1] - 2.0*f[j] + f[j-1])/(dmu*dmu)

    dfdmumu[0]  = (f[2] - 2.0*f[1] + f[0])/(dmu*dmu)
    dfdmumu[-1] = (f[-3] - 2.0*f[-2] + f[-1])/(dmu*dmu)
    
    for j in range(0,len(mu)):
        dfdt[j]    = Dmumu * ( -2*mu[j]*dfdmu[j] + (1.0 - mu[j]*mu[j])*dfdmumu[j])
        f[j] = f[j] + dfdt[j]*dt

    fmin = np.amin(f)
    fmax = np.amax(f)

    A[t] = 1 - np.log(fmax/fmin)/(v**2)
      
    fplot =  np.exp(- vpara*vpara - A[t] * vperp*vperp)

    #if time[t] == tplot[int(t/10)]:

    #    ax_fmu = plt.subplot(221)
    #    ax_fmu.plot(mu,f,color='black')
    #    ax_fmu.set_xlim(-1.0,1.0)
    #    ax_fmu.set_ylim(0.0,0.001)
    #    ax_fmu.set_xlabel('$\\mu$',fontsize=15)
    #    ax_fmu.set_ylabel('$f$',fontsize=15)
    #    
    #    ax_fv = plt.subplot(222)
    #    ax_fv.pcolormesh(vpara,vperp,fplot,cmap='viridis',shading='auto')
    #    ax_fv.set_xlim(vrange[0]/3,vrange[-1]/3)
    #    ax_fv.set_ylim(vrange[0]/3,vrange[-1]/3)
    #    ax_fv.set_xlabel('$v_\\parallel$',fontsize=15)
    #    ax_fv.set_ylabel('$v_\\bot$',fontsize=15)
    #    ax_fv.set_aspect('equal')
    #    
    #    ax_A = plt.subplot(212)
    #    ax_A.plot(time[:t],A[:t],color='black')
    #    ax_A.set_xlim(0.0,12.0)
    #    ax_A.set_ylim(A[0],1.0)
    #    ax_A.set_xlabel('Time (s)',fontsize=15)
    #    ax_A.set_ylabel('A',fontsize=15)
    #    
    #    plt.subplots_adjust(hspace = 0.4)
    #    
    #    plt.savefig(path_fig+'A_v'+str(v)+'_D'+str(Dmumu)+'_t'+"{:.3f}".format(time[t])+'.png',,dpi=300)
    #    print(path_fig+'A_v'+str(v)+'_D'+str(Dmumu)+'_t'+"{:.3f}".format(time[t])+'.png')
    #    plt.close()

print(A[-1])
print(f)

np.save(path_save+'A_'+str(Ainit)+'_v'+str(v)+'_D'+str(Dmumu)+'_dt'+str(dt)+'.npy',A)
print(path_save+'A_'+str(Ainit)+'_v'+str(v)+'_D'+str(Dmumu)+'_dt'+str(dt)+'.npy')

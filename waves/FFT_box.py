import pytools as pt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
import scipy
import scipy.ndimage
import matplotlib.colors as colors
import imageio
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

run = sys.argv[1]

bulkStart = int(sys.argv[2])
bulkEnd   = int(sys.argv[3])
timetot   = bulkEnd - bulkStart + 1

fontsize = 20

if run == '300':
    path_bulk = '/wrk-vakka/users/dubart/diff_test/dmumu/'+run+'/bulk/'
elif run == 'Markus':
    path_bulk = '/wrk-vakka/users/markusb/weirddiffusion/Losscone_may31/'
else:
    path_bulk = '/wrk-vakka/users/dubart/diff_test/dmumu/diff/'+run+'/CFL1.0/bulk/'

path_fig = '/wrk-vakka/users/dubart/analysis/waves/fig/'

x_length = 32
y_length = 32
z_length = 1
z_min = 0
z_max = 3e5
x_min = 0.0
x_max = x_length*z_max
y_min = 0.0
y_max = y_length*z_max

x_coordinates_para = np.linspace(x_min+z_max,x_max-z_max,x_length)
y_coordinates_para = y_max / 2.0
z_coordinates_para = z_max / 2.0

x_coordinates_perp = x_max / 2.0
y_coordinates_perp = np.linspace(y_min+z_max,y_max-z_max,y_length)
z_coordinates_perp = z_max / 2.0

B_kpara = np.zeros([len(x_coordinates_para),3,timetot])
#E_kpara = np.zeros([len(x_coordinates_para),3,timetot])
B_kperp = np.zeros([len(y_coordinates_perp),3,timetot])
#E_kperp = np.zeros([len(y_coordinates_perp),3,timetot])
v_kpara = np.zeros([len(x_coordinates_para),3,timetot])
v_kperp = np.zeros([len(y_coordinates_perp),3,timetot])
n_kpara = np.zeros([len(y_coordinates_perp),timetot])
n_kperp = np.zeros([len(y_coordinates_perp),timetot])

kpara_omega_B = np.zeros([timetot,3,len(x_coordinates_para)])
kperp_omega_B = np.zeros([timetot,3,len(y_coordinates_perp)])

time = np.zeros(timetot)

#Constants (SI units)
q        = 1.60217e-19
mu0      = 1.256637061e-6
c        = 299792458
mp       = 1.6726219e-27
kB       = 1.38064852e-23
me       = 9.1093826e-31
epsilon0 = 8.85418781762e-12
gamma    = 5.0/3.0


for step in range(bulkStart,bulkEnd+1):

    if step == 0:
        bulkname = 'initial-grid.0000000.vlsv'
    else:
        bulkname = 'bulk.'+str(step).rjust(7,'0')+'.vlsv'

    f = pt.vlsvfile.VlsvReader(path_bulk + bulkname)

    time[step] = f.read_parameter('time')

    xsize = f.read_parameter('xcells_ini')
    xmax  = f.read_parameter("xmax")
    xmin  = f.read_parameter("xmin")
    dx    = (xmax-xmin)/xsize

    for i in range(0,len(x_coordinates_para)):
        cellid_para        = int(f.get_cellid([x_coordinates_para[i],y_coordinates_para,z_coordinates_para]))
        B_kpara[i,:,step]  = f.read_variable('vg_b_vol',cellid_para)
        #E_kpara[i,:,step]  = f.read_variable('vg_e_vol',cellid_para)
        
        v_kpara[i,:,step] = f.read_variable('vg_v',cellid_para)

        n_kpara[i,step] = f.read_variable('vg_rho',cellid_para)

    for i in range(0,len(y_coordinates_perp)):
        cellid_perp        = int(f.get_cellid([x_coordinates_perp,y_coordinates_perp[i],z_coordinates_perp]))
        B_kperp[i,:,step]  = f.read_variable('vg_b_vol',cellid_perp)
        #E_kperp[i,:,step]  = f.read_variable('vg_e_vol',cellid_perp)

        v_kperp[i,:,step] = f.read_variable('vg_v',cellid_perp)

        n_kperp[i,step] = f.read_variable('vg_rho',cellid_perp)
    
B_kpara_norm = np.linalg.norm(np.mean(np.mean(B_kpara,axis=2),axis=0)) 
B_kperp_norm = np.linalg.norm(np.mean(np.mean(B_kperp,axis=2),axis=0)) 
v_kpara_norm = np.linalg.norm(np.mean(np.mean(v_kpara,axis=2),axis=0)) 
v_kperp_norm = np.linalg.norm(np.mean(np.mean(v_kperp,axis=2),axis=0)) 
n_kpara_avg  = np.mean(np.mean(n_kpara,axis=1))
n_kperp_avg  = np.mean(np.mean(n_kperp,axis=1))

d_i_kpara  = np.sqrt(mp*epsilon0*c*c/(n_kpara_avg*q*q))
d_i_kperp  = np.sqrt(mp*epsilon0*c*c/(n_kperp_avg*q*q))
w_ci_kpara = q*B_kpara_norm/mp
w_ci_kperp = q*B_kperp_norm/mp
vA         = B_kpara_norm / np.sqrt(mu0 * n_kpara_avg * mp)

vpara = np.dot(np.mean(np.mean(B_kpara,axis=2),axis=0),np.mean(np.mean(v_kpara,axis=2),axis=0)) / B_kpara_norm

for i in range(0,3):

    kpara_omega_B[:,i,:] = np.fft.fftshift(np.fft.fft2(B_kpara[:,i,:].T))
    kperp_omega_B[:,i,:] = np.fft.fftshift(np.fft.fft2(B_kperp[:,i,:].T))
 
dt = abs(time[0] - time[1])

omegamin = 0.0
omegamax = np.pi / dt / ((w_ci_kpara + w_ci_kperp)/2)

kmin = -np.pi / dx / ((d_i_kpara + d_i_kperp)/2)
kmax = np.pi / dx / ((d_i_kpara + d_i_kperp)/2)

krange = np.linspace(kmin, kmax, 2000)

awaves_para = abs(vA * (abs(krange)/d_i_kpara) + vpara * (krange/d_i_kpara)) / w_ci_kpara

# Plot
fig, (ax_para, ax_perp)    = plt.subplots(1,2,sharey=True)

image = ax_para.imshow(abs(kpara_omega_B[int((kpara_omega_B.shape[0])/2):,1,:]), norm=colors.LogNorm(), extent = [kmin,kmax,omegamin,omegamax], aspect = 'auto')

ax_para.set_xlim(kmin, kmax)
ax_para.set_ylim(omegamin, omegamax)
ax_para.set_xlabel("$k_\\parallel * d_i$",fontsize=fontsize)
ax_para.set_ylabel("$\\frac{\\omega}{\\Omega_c}$",fontsize=fontsize)

ax_para.plot(krange,awaves_para,linewidth=2,color='blue',label='Alfven Waves')

image = ax_perp.imshow(abs(kperp_omega_B[int((kperp_omega_B.shape[0])/2):,1,:]), norm=colors.LogNorm(), extent = [kmin,kmax,omegamin,omegamax], aspect = 'auto')

ax_perp.set_xlim(kmin, kmax)
ax_perp.set_ylim(omegamin, omegamax)
ax_perp.set_xlabel("$k_\\bot * d_i$",fontsize=fontsize)
ax_perp.set_ylabel("$\\frac{\\omega}{\\Omega_c}$",fontsize=fontsize)

plt.savefig(path_fig+'FFT_'+run+'.png',dpi=300)
print(path_fig+'FFT_'+run+'.png')

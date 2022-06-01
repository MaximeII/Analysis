import pytools as pt
import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from multiprocessing import Pool
from matplotlib import rc, rcParams
from matplotlib.ticker import LogLocator


# Avoids opening a figure window
if str(matplotlib.get_backend()) != 'Agg':
    plt.switch_backend('Agg') 

numproc = 20

run = sys.argv[1]

path_bulk = '/wrk/users/dubart/diff_test/proc_test/'+run+'/'
path_fig  = '/wrk/users/dubart/analysis/subgrid/300_test/Green_diff/fig/proc_test/test_files/'+run+'/'
path_data = '/wrk/users/dubart/analysis/subgrid/300_test/Green_diff/data/'

if run == 'mu_bulk':
    path_save = '/wrk/users/dubart/diff_test/proc_test/mu_files/'
    figname   = 'mu_files'
    figtitle  = '$\\mu$ diffusion'
    varmin    = -1.0
    varmax    = 1.0
else:
    path_save = '/wrk/users/dubart/diff_test/proc_test/cart_files/'
    figname   = 'cart_files'
    figtitle  = 'Cartesian diffusion'
    varmin    = -1.0
    varmax    = 1.0

RE = 6371e3 # m

vmin = -2.0
vmax = 2.0

labelsize = 20
ticksize = 20

Ddtmin = 1e-6
Ddtmax = 0.01
tmin   = 0.0
tmax   = 15.0

# Figure layout
mh = 0.03 
mv = 0.15
ecb = 0.015
lcb = 0.02
ev = 0.1
h1 = 1-1.5*mv
h2 = (h1-ev)/2.
l1 = 0.32
l2 = l1*h2/h1
eh = (1.-l1-2*l2-2.5*ecb-2*lcb-2*mh)/3.
fontscale = 2.5
figratio = h2/l2
deltaX = mh+l1+ecb+lcb+2*eh-0.05


if len(sys.argv) == 4: # Starting and end frames given
    timetot = range(int(sys.argv[2]), int(sys.argv[3]),1)
else: # Only starting frame given, generate one frame
    timetot = range(int(sys.argv[2]), int(sys.argv[2])+1, 1)
for step in timetot:

    fig = plt.figure(figsize=(8.*figratio,8.))

    if step == 0:
        bulkname       = "initial-grid.0000000.vlsv"
        if run == 'mu_bulk':
            panel_fmu      = '(a)'
            panel_dfdtmu    = '(b)'
            panel_dfdt_per = '(c)'
            panel_dfdt_par = '(d)'
            panel_f_per    = '(e)'
            panel_f_par    = '(f)'
        else:
            panel_dfdt_per = '(a)'
            panel_dfdt_par = '(b)'
            panel_f_per    = '(c)'
            panel_f_par    = '(d)'
    else:
        bulkname = "bulk."+str(step).rjust(7,'0')+".vlsv"
        if run == 'mu_bulk':
            panel_fmu      = '(g)'
            panel_dfdtmu    = '(h)'
            panel_dfdt_per = '(i)'
            panel_dfdt_par = '(j)'
            panel_f_per    = '(k)'
            panel_f_par    = '(l)'
        else:
            panel_dfdt_per = '(e)'
            panel_dfdt_par = '(f)'
            panel_f_per    = '(g)'
            panel_f_par    = '(h)'


    cid = 1

    f = pt.vlsvfile.VlsvReader(path_bulk+bulkname)
    coords = f.get_cell_coordinates(cid)

    # VDFs and their colorbar -----------------------------------------------
    ax_VDF_per = fig.add_axes([deltaX+l2+eh,mv+h2+ev,l2,h2])
    ax_VDF_par = fig.add_axes([deltaX+l2+eh,mv-0.04,l2,h2])
    cax_VDF    = fig.add_axes([deltaX+2*l2+eh+1.5*ecb,mv-0.04,lcb,h1+0.04])

    pt.plot.plot_vdf(filename   = path_bulk+bulkname,
                     cellids    = [cid],
                     xy         = 1,
                     box        = [-2.0e6,2.0e6,-2.0e6,2.0e6],
                     fmin       = 1e-15, fmax = 1e-11,
                     cbulk      = 0,
                     axes       = ax_VDF_par,cbaxes=cax_VDF,
                     slicethick = 1,
                     axisunit   = 6,
                     title      = '',
                     scale      = 2.5,
                     cbtitle    = '',
                     #setThreshold = 1e-21,
                     colormap   = 'viridis')

    pt.plot.plot_vdf(filename   = path_bulk+bulkname,
                     cellids    = [cid],
                     yz         = 1,
                     box        = [-2.0e6,2.0e6,-2.0e6,2.0e6],
                     fmin       = 1e-15, fmax = 1e-11,
                     cbulk      = 0,
                     axes       = ax_VDF_per,nocb=1,
                     slicethick = 1,
                     axisunit   = 6,
                     title      = '',
                     scale      = 2.5,
                     #setThreshold = 1e-21,
                     colormap   = 'viridis')
   
    cb_text = r'$f($v$)$ [m$^{-6}$ s$^{3}$]'
    plt.text(-0.1,-0.08,cb_text,transform=cax_VDF.transAxes,fontsize=20,weight='bold')
    ax_VDF_per.text(-1.9,1.6,panel_f_per,fontsize=22,weight='bold')
    ax_VDF_par.text(-1.9,1.6,panel_f_par,fontsize=22,weight='bold')
   

    # dfdt axes -------------------------------------------------------------
    ax_dfdt_par = fig.add_axes([deltaX-1.5*ecb-lcb-0.05,mv-0.04,l2,h2])
    ax_dfdt_per = fig.add_axes([deltaX-1.5*ecb-lcb-0.05,mv+h2+ev,l2,h2])
    dfdt_cax    = fig.add_axes([deltaX+l2-0.03-0.04,mv-0.04,lcb,h1+0.04]) 
 
    # Get dfdt data
    data  = np.loadtxt(path_save+'dfdt_array_'+str(step).rjust(7,'0')+'.txt')
    VX    = data[:,0]
    VY    = data[:,1]
    VZ    = data[:,2]
    dfdt  = data[:,3]
    fcell = data[:,4]
   
    ratio = dfdt/fcell

    minX = min(VX)
    maxX = max(VX)
    minY = min(VY)
    maxY = max(VY)
    minZ = min(VZ)
    maxZ = max(VZ)

    # Build back 3D array of VDF
    x = np.arange(minX,maxX+30000,30000)
    y = np.arange(minY,maxY+30000,30000)
    z = np.arange(minZ,maxZ+30000,30000)

    var3D = np.zeros([len(x),len(y),len(z)])

    for i in range(0,len(dfdt)):
        indX = np.where(x == VX[i])[0][0]
        indY = np.where(y == VY[i])[0][0]
        indZ = np.where(z == VZ[i])[0][0]
        #var3D[indX,indY,indZ] = dfdt[i]
        var3D[indX,indY,indZ] = ratio[i]

    # Build slice at VX,VY,VZ = 0 (avg of +-)    
    varXY = np.mean(var3D[:,:,np.where(z == -15000)[0][0]:np.where(z == 15000)[0][0]+1],axis=2)
    
    varXZ = np.mean(var3D[:,np.where(y == -15000)[0][0]:np.where(y == 15000)[0][0]+1,:],axis=1)
    
    varYZ = np.mean(var3D[np.where(x == -15000)[0][0]:np.where(x == 15000)[0][0]+1,:,:],axis=0)
    
    XY, YX = np.meshgrid(x/1e6,y/1e6)
    XZ, ZX = np.meshgrid(x/1e6,z/1e6)
    YZ, ZY = np.meshgrid(y/1e6,z/1e6)
    
    # Plot dfdt
    figdfdt = ax_dfdt_per.pcolormesh(YZ,ZY,varYZ.T,cmap='seismic',shading='auto',vmin=varmin,vmax=varmax)
    ax_dfdt_par.pcolormesh(XY,YX,varXY.T,cmap='seismic',shading='auto',vmin=varmin,vmax=varmax)

    ax_dfdt_par.set_xlim(vmin,vmax)
    ax_dfdt_par.set_ylim(vmin,vmax)
    ax_dfdt_par.set_xlabel(r'$v_x$ [$10^6$ m s$^{-1}$]',fontsize=labelsize,weight='bold')
    ax_dfdt_par.set_ylabel(r'$v_y$ [$10^6$ m s$^{-1}$]',fontsize=labelsize,weight='bold') 
    ax_dfdt_par.xaxis.set_tick_params(labelsize=ticksize)
    ax_dfdt_par.yaxis.set_tick_params(labelsize=ticksize)
    
    ax_dfdt_per.set_xlim(vmin,vmax)
    ax_dfdt_per.set_ylim(vmin,vmax)
    ax_dfdt_per.set_xlabel(r'$v_y$ [$10^6$ m s$^{-1}$]',fontsize=labelsize,weight='bold')
    ax_dfdt_per.set_ylabel(r'$v_z$ [$10^6$ m s$^{-1}$]',fontsize=labelsize,weight='bold') 
    ax_dfdt_per.xaxis.set_tick_params(labelsize=ticksize)
    ax_dfdt_per.yaxis.set_tick_params(labelsize=ticksize)

    cb_dfdt = plt.colorbar(figdfdt,cax=dfdt_cax)
    cb_dfdt.ax.tick_params(labelsize=20)
    plt.text(-0.05,-0.08,'$\\partial_t f($v$)$/$f($v$)$',transform=dfdt_cax.transAxes,fontsize = 23,weight='bold') 
    ax_dfdt_per.text(-1.9,1.6,panel_dfdt_per,fontsize=22,weight='bold')
    ax_dfdt_par.text(-1.9,1.6,panel_dfdt_par,fontsize=22,weight='bold')
    ax_dfdt_per.grid('grey',linestyle='-',lw=1.0)
    ax_dfdt_par.grid('grey',linestyle='-',lw=1.0)
    ax_dfdt_per.tick_params(axis='x',which='minor')
    ax_dfdt_par.tick_params(axis='x',which='minor')
    ax_dfdt_per.tick_params(axis='y',which='minor')
    ax_dfdt_par.tick_params(axis='y',which='minor')
    ax_dfdt_par.xaxis.set_ticks(np.arange(vmin,vmax+1.0,1.0))
    ax_dfdt_par.yaxis.set_ticks(np.arange(vmin,vmax+0.5,0.5))
    ax_dfdt_per.xaxis.set_ticks(np.arange(vmin,vmax+1.0,1.0))
    ax_dfdt_per.yaxis.set_ticks(np.arange(vmin,vmax+0.5,0.5))

    labels = ax_dfdt_par.get_xticklabels() + ax_dfdt_par.get_yticklabels() + ax_dfdt_per.get_xticklabels() + ax_dfdt_per.get_yticklabels()
    [label.set_fontweight('bold') for label in labels]

    if run == 'mu_bulk': 

        minmu = -1.0
        maxmu = 1.0
        minv  = 30000.0/1e6
        maxv  = 1.39257e+07/1e6

        # dfdtmu axes ---------------------------------------------------------
        ax_dfdtmu    = fig.add_axes([deltaX-2.0*(1.5*ecb+lcb)-l2-eh-0.11,mv-0.04,l2,h2])
        dfdtmu_cax   = fig.add_axes([deltaX-1.5*ecb-lcb-eh-0.03-0.1,mv-0.04,lcb,h2])

        # Get dftdmu data
        data = np.loadtxt(path_save+'dfdt_mu_array_'+str(step).rjust(7,'0')+'.txt')
    
        murange = np.linspace(minmu,maxmu,data.shape[1]+1)
        vrange  = np.linspace(minv,maxv,data.shape[0]+1)

        muv, vmu = np.meshgrid(murange, vrange)

        dfdtmumin = -0.2
        dfdtmumax = 0.2

        # Plot dfdtmu
        figdfdtmu = ax_dfdtmu.pcolormesh(muv, vmu, data,vmin=dfdtmumin,vmax=dfdtmumax,cmap='seismic',shading='auto')
        ax_dfdtmu.set_xlim(minmu,maxmu)
        ax_dfdtmu.set_ylim(0.0,2.0)
        ax_dfdtmu.set_xlabel('$\\mu$',fontsize=labelsize)
        ax_dfdtmu.set_ylabel(r'$V$ [$10^6$ m s$^{-1}$]',fontsize=labelsize,weight='bold')   
        ax_dfdtmu.xaxis.set_tick_params(labelsize=ticksize)
        ax_dfdtmu.yaxis.set_tick_params(labelsize=ticksize)

        ax_dfdtmu.xaxis.set_ticks(np.arange(minmu,maxmu+0.5,0.5))
        ax_dfdtmu.yaxis.set_ticks(np.arange(0.0,2.5,0.5)) 

        cb_dfdtmu = plt.colorbar(figdfdtmu,cax=dfdtmu_cax)
        cb_dfdtmu.ax.tick_params(labelsize=20)
        plt.text(-0.06,-0.18,'$\\partial_t f(\\mu,v)$ [m$^{-4}$]',transform=dfdtmu_cax.transAxes,fontsize = 23,weight='bold')    
        ax_dfdtmu.text(-0.95,1.8,panel_dfdtmu,fontsize=22,weight='bold')
 
        # fmu axes -----------------------------------------------------------
        ax_fmu  = fig.add_axes([deltaX-2.0*(1.5*ecb+lcb)-l2-eh-0.11,mv+h2+ev,l2,h2])
        fmu_cax = fig.add_axes([deltaX-1.5*ecb-lcb-eh-0.03-0.1,mv+h2+ev,lcb,h2])

        # Get fmu data
        data = np.loadtxt(path_save+'muv_array_'+str(step).rjust(7,'0')+'.txt')
     
        fmumin = 0.0
        fmumax = 4.0
           
        # Plot fmu
        figfmu = ax_fmu.pcolormesh(muv, vmu, data,vmin=fmumin,vmax=fmumax,cmap='viridis',shading='auto')
        ax_fmu.set_xlim(minmu,maxmu)
        ax_fmu.set_ylim(0.0,2.0)
        ax_fmu.set_xlabel('$\\mu$',fontsize=labelsize)
        ax_fmu.set_ylabel(r'$V$ [$10^6$ m s$^{-1}$]',fontsize=labelsize,weight='bold')
        ax_fmu.xaxis.set_tick_params(labelsize=ticksize)
        ax_fmu.yaxis.set_tick_params(labelsize=ticksize)

        ax_fmu.xaxis.set_ticks(np.arange(minmu,maxmu+0.5,0.5))
        ax_fmu.yaxis.set_ticks(np.arange(0.0,2.5,0.5)) 

        cb_fmu = plt.colorbar(figfmu,cax=fmu_cax)
        cb_fmu.ax.tick_params(labelsize=20)

        plt.text(-0.06,-0.18,'$f(\\mu,v)$ [m$^{-4}$ s$^{1}$]',transform=fmu_cax.transAxes,fontsize = 23,weight='bold')     
        ax_fmu.text(-0.95,1.8,panel_fmu,fontsize=22,weight='bold',color='white')

        labels = ax_dfdtmu.get_xticklabels() + ax_dfdtmu.get_yticklabels() + ax_fmu.get_xticklabels() + ax_fmu.get_yticklabels()
        [label.set_fontweight('bold') for label in labels]

    #plt.figtext(0.02,0.5,figtitle+', t = '+str(step/100.0)+'s',fontsize=20,weight='bold') 
    plt.suptitle(figtitle+', t = '+str(step/100.0)+'s',fontsize=20,weight='bold')
    plt.savefig(path_fig+figname+'_'+str(step).rjust(7,'0')+'.png')
    print(path_fig+figname+'_'+str(step).rjust(7,'0')+'.png')

    plt.close()


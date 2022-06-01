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
xmargin  = 0.12
ymargin  = 0.11
axwidth  = 0.25
axheight = 0.35
cbmargin = 0.015
cbwidth  = 0.02
axspace  = 0.06

if len(sys.argv) == 4: # Starting and end frames given
    timetot = range(int(sys.argv[2]), int(sys.argv[3]),1)
else: # Only starting frame given, generate one frame
    timetot = range(int(sys.argv[2]), int(sys.argv[2])+1, 1)
for step in timetot:

    fig = plt.figure(figsize=(12.,8.))

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
    ax_VDF_per = fig.add_axes([xmargin+axwidth+cbmargin+cbwidth+axspace*3+0.01,ymargin+axheight+axspace*2,axwidth,axheight])
    ax_VDF_par = fig.add_axes([xmargin+axwidth+cbmargin+cbwidth+axspace*3+0.01,ymargin,axwidth,axheight])
    cax_VDF    = fig.add_axes([xmargin+axwidth+cbmargin+cbwidth+axspace*3+axwidth+cbmargin+0.01,ymargin,cbwidth,axheight*2+axspace*2])

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
    ax_VDF_par.text(2.0,-3.0,cb_text,fontsize=20,weight='bold')
    ax_VDF_per.text(-1.9,1.6,panel_f_per,fontsize=22,weight='bold')
    ax_VDF_par.text(-1.9,1.6,panel_f_par,fontsize=22,weight='bold')
   

    # dfdt axes -------------------------------------------------------------
    ax_dfdt_per = fig.add_axes([xmargin,ymargin+axheight+axspace*2,axwidth,axheight])
    ax_dfdt_par = fig.add_axes([xmargin,ymargin,axwidth,axheight])
    dfdt_cax    = fig.add_axes([xmargin+axwidth+cbmargin,ymargin,cbwidth,axheight*2+axspace*2])

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
    plt.text(-0.07,-0.11,'$\\partial_t f($v$)$/$f($v$)$',transform=dfdt_cax.transAxes,fontsize = 23,weight='bold') 
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

    #plt.figtext(0.02,0.5,figtitle+', t = '+str(step/100.0)+'s',fontsize=20,weight='bold') 
    plt.suptitle(figtitle+', t = '+str(step/100.0)+'s',fontsize=20,weight='bold')
    plt.savefig(path_fig+figname+'_'+str(step).rjust(7,'0')+'.png')
    print(path_fig+figname+'_'+str(step).rjust(7,'0')+'.png')

    plt.close()


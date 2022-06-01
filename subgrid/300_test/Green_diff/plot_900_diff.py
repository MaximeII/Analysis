import pytools as pt
import sys, os, socket
import numpy as np
import matplotlib.pyplot as plt
import colormaps as cmaps
import matplotlib
import matplotlib.patches as patches


# Avoids opening a figure window
if str(matplotlib.get_backend()) != 'Agg':
    plt.switch_backend('Agg')


fileLocation='/wrk/users/dubart/diff_test/900km/diffusion/restart/bulk/'

panel_VDF_par1 = '(a)'
panel_VDF_par2 = '(b)'
panel_VDF_perp = '(c)'

scale_vdf = 2.5

outputLocation=outputdir='/wrk/users/dubart/analysis/subgrid/300_test/Green_diff/fig/proc_test/Colormap/'
path_save = '/wrk/users/dubart/analysis/subgrid/300_test/Green_diff/data/'


if len(sys.argv) == 3: # Starting and end frames given
    timetot = range(int(sys.argv[1]), int(sys.argv[2]),1)
else: # Only starting frame given, generate one frame
    timetot = range(int(sys.argv[1]), int(sys.argv[1])+1, 1)
for step in timetot:

    if step == 0:
        bulkname = "initial-grid.0000000.vlsv"       
    else:
        # Source data file
        bulkname = "bulk."+str(step).rjust(7,'0')+".vlsv"
        #bulkname = "restart."+str(step).rjust(7,'0')+".2021-11-26_18-31-02.vlsv"
  
    if step == 837:
        panel_VDF_par1 = '(a)'
        panel_VDF_par2 = '(b)'
        panel_VDF_perp = '(c)'
    else:
        panel_VDF_par1 = '(d)'
        panel_VDF_par2 = '(e)'
        panel_VDF_perp = '(f)'

    print(bulkname)

    # CellIDs of interest
    cid1 = int(np.load(path_save+'CellID_diffusion.npy'))
    Re = 6371e3 # m

    f = pt.vlsvfile.VlsvReader(fileLocation+bulkname)
    coords = f.get_cell_coordinates(cid1)

    # Figure Layout (in pxl from 0-1)
    xmargin  = 0.05
    ymargin  = 0.12
    axwidth  = 0.25
    axheight = 0.8
    cbmargin = 0.015
    cbwidth  = 0.02
    axspace  = 0.03
    figratio = 2.421875

    fig = plt.figure(figsize=(22,6.5))

    ax_VDF_par1 = fig.add_axes([xmargin,ymargin,axwidth,axheight])
    ax_VDF_par2 = fig.add_axes([xmargin+axwidth+axspace*2,ymargin,axwidth,axheight])
    ax_VDF_perp = fig.add_axes([xmargin+axwidth+axspace*2+axwidth+axspace*2,ymargin,axwidth,axheight])
    ax_cb       = fig.add_axes([xmargin+axwidth+axspace*2+axwidth+axspace*2+axwidth+cbmargin,ymargin,cbwidth,axheight])

    VDF_min = -2.0e6
    VDF_max = 2.0e6

    pt.plot.plot_vdf(filename=fileLocation+bulkname,
                     cellids=[cid1],
                     xy=1,
                     box=[VDF_min,VDF_max,VDF_min,VDF_max],
                     fmin=1e-15, fmax=1e-11,
                     slicethick=1,
                     axes=ax_VDF_par1,nocb=1,
                     cbulk=1,
                     axisunit=6,
                     title='',
                     cbtitle='',
                     scale=scale_vdf,
                     colormap='viridis')

    pt.plot.plot_vdf(filename=fileLocation+bulkname,
                     cellids=[cid1],
                     xz=1,
                     box=[VDF_min,VDF_max,VDF_min,VDF_max],
                     fmin=1e-15, fmax=1e-11,
                     slicethick=1,
                     axes=ax_VDF_par2,nocb=1,
                     cbulk=1,
                     axisunit=6,
                     title='',
                     cbtitle='',
                     scale=scale_vdf,
                     colormap='viridis')

    pt.plot.plot_vdf(filename=fileLocation+bulkname,
                     cellids=[cid1],
                     yz=1,
                     box=[VDF_min,VDF_max,VDF_min,VDF_max],
                     fmin=1e-15, fmax=1e-11,
                     slicethick=1,
                     axes=ax_VDF_perp,cbaxes=ax_cb,
                     cbulk=1,
                     axisunit=6,
                     title='',
                     cbtitle='',
                     scale=scale_vdf,
                     colormap='viridis')

    plt.suptitle('t = '+str(step/2.0)+'s',fontsize=20,weight='bold')
    plt.text(-1.0,-0.12,'$f($v$)$ [m$^{-6}$ s$^{3}$]',transform=ax_cb.transAxes,fontsize=20,weight='bold')
    ax_VDF_par1.text(-1.95,1.8,panel_VDF_par1,fontsize=22,weight='bold')
    ax_VDF_par2.text(-1.95,1.8,panel_VDF_par2,fontsize=22,weight='bold')
    ax_VDF_perp.text(-1.95,1.8,panel_VDF_perp,fontsize=22,weight='bold')

    plt.savefig(outputLocation+"900_diff_"+str(step).rjust(7,'0')+".png",dpi=300)
    print(outputLocation+"900_diff_"+str(step).rjust(7,'0')+".png")
    plt.close()

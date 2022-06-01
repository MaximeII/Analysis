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

run       = sys.argv[1]

scale_map = 2.5
scale_vdf = 2.5

RE = 6371e3

#fileLocation='/wrk/users/dubart/300_test/'+run+'/mu_bulk/'
#fileLocation='/wrk/users/dubart/300_test/'+run+'/restart/'
#fileLocation = '/wrk/users/dubart/300_test/proc_test/'+run+'/'
#fileLocation='/wrk/users/dubart/300_test/900km/diffusion/restart/bulk/'
fileLocation='/wrk/users/dubart/diff_test/proc_test/mu_bulk/'

panel = '(a)'
panel_VDF1 = '(d)'
panel_VDF2 = '(e)'
panel_VDF3 = '(f)'

#outputLocation=outputdir='/wrk/users/dubart/analysis/subgrid/300_test/Green_diff/fig/'+run+'/restart/'
outputLocation=outputdir='/wrk/users/dubart/analysis/subgrid/300_test/Green_diff/fig/proc_test/Colormap/'
path_save = '/wrk/users/dubart/analysis/subgrid/300_test/Green_diff/data/'


if len(sys.argv) == 4: # Starting and end frames given
    timetot = range(int(sys.argv[2]), int(sys.argv[3]),1)
else: # Only starting frame given, generate one frame
    timetot = range(int(sys.argv[2]), int(sys.argv[2])+1, 1)
for step in timetot:

    if step == 0:
        bulkname = "initial-grid.0000000.vlsv"
    else:
        # Source data file
        bulkname = "bulk."+str(step).rjust(7,'0')+".vlsv"
        #bulkname = "restart."+str(step).rjust(7,'0')+".2021-11-26_18-31-02.vlsv"
    print(bulkname)

    # CellIDs of interest
    #cid1 = int(np.load(path_save+'CellID_'+run+'.npy'))
    cid1 = 1
    Re = 6371e3 # m

    f = pt.vlsvfile.VlsvReader(fileLocation+bulkname)
    coords = f.get_cell_coordinates(cid1)

    # Figure layout
    mh = 0.05
    mv = 0.12
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

    fig = plt.figure(figsize=(8.*figratio,8.))
    # Colormap and its colorbar
    ax_col = fig.add_axes([mh,mv,l1,h1])
    cax_col = fig.add_axes([mh+l1+ecb,mv,lcb,h1])
    # VDFs and their colorbar
    deltaX = mh+l1+ecb+lcb+2*eh
    ax_VDF1_par  = fig.add_axes([deltaX,mv+h2+ev,l2,h2])
    ax_VDF1_per  = fig.add_axes([deltaX+l2+eh,mv+h2+ev,l2,h2])
    ax_VDF1_per2 = fig.add_axes([deltaX,mv,l2,h2])
   # ax_VDF2_per = fig.add_axes([deltaX+l2+eh,mv,l2,h2])
    cax_VDF = fig.add_axes([deltaX+2*l2+eh+1.5*ecb,mv,lcb,h1-0.01])


    # Colormap with fg_b
    # pt.plot.plot_colormap(filename=fileLocation+bulkname,
    #                       #var='vg_b_vol',
    #                       #vmin=1.0e-8,vmax=2.0e-8,
    #                       var='proton/vg_rho',
    #                       vmin=1.0e6,vmax=5.0e6,
    #                       boxm=[coords[0]-9e5,coords[0]+9e5,coords[1]-9e5,coords[1]+9e5],
    #                       run=run,
    #                       colormap='seismic',
    #                       scale=scale_map,
    #                       axes=ax_col,cbaxes=cax_col,
    #                       cbtitle = '',
    #                       outputdir=outputLocation,lin=1)

    # Addition of location of cells with VDFs
    vlsvReader=pt.vlsvfile.VlsvReader(fileLocation+bulkname)
    xCid,yCid,zCid = vlsvReader.get_cell_coordinates(cid1)
    VDFcoord_1 = [xCid/Re,yCid/Re,zCid/Re]

    ax_col.scatter(VDFcoord_1[0], VDFcoord_1[1], color='black',marker='o',s=150)
    ax_col.scatter(VDFcoord_1[0], VDFcoord_1[1], color='white',marker='o',s=70)

    VDF_min = -2.0e6
    VDF_max = 2.0e6

    # VDFs of first cell
    pt.plot.plot_vdf(filename=fileLocation+bulkname,
                     cellids=[cid1],
                     xy=1,
                     box=[VDF_min,VDF_max,VDF_min,VDF_max],
                     fmin=1e-15, fmax=1e-11,
                     slicethick=1,
                     axes=ax_VDF1_par,cbaxes=cax_VDF,
                     cbulk=1,
                     #center = [0.0,1.0e6,0.0],
                     axisunit=6,
                     title='',
                     scale=scale_vdf,
                     cbtitle='',
                     colormap='viridis')

    pt.plot.plot_vdf(filename=fileLocation+bulkname,
                     cellids=[cid1],
                     yz=1,
                     box=[VDF_min,VDF_max,VDF_min,VDF_max],
                     fmin=1e-15, fmax=1e-11,
                     slicethick=1,
                     axes=ax_VDF1_per,nocb=1,
                     cbulk=1,
                     #center = [0.0,1.0e6,0.0],
                     axisunit=6,
                     title='',
                     scale=scale_vdf,
                     colormap='viridis')

    pt.plot.plot_vdf(filename=fileLocation+bulkname,
                     cellids=[cid1],
                     xz=1,
                     box=[VDF_min,VDF_max,VDF_min,VDF_max],
                     fmin=1e-15, fmax=1e-11,
                     slicethick=1,
                     axes=ax_VDF1_per2,nocb=1,
                     cbulk=1,
                     #center = [0.0,1.0e6,0.0],
                     axisunit=6,
                     title='',
                     scale=scale_vdf,
                     colormap='viridis')

    plt.text(0.28, 1.1,panel,transform=ax_col.transAxes,fontsize=22,weight='bold')
    plt.text(-1.0,-0.08,'$f(v)$ [m$^{-6}$ s$^{3}$]',transform=cax_VDF.transAxes,fontsize=20,weight='bold')
    #plt.text(-0.1,-0.08,'$n_\\mathrm{proton}$ [m$^{-3}$]',transform=cax_col.transAxes,fontsize=20,weight='bold')
    ax_col.text(-1.3,21.1,panel,fontsize=25,weight='bold')
    ax_VDF1_par.text(-1.9,1.6,panel_VDF1,fontsize=22,weight='bold')
    ax_VDF1_per.text(-1.9,1.6,panel_VDF2,fontsize=22,weight='bold')
    ax_VDF1_per2.text(-1.9,1.6,panel_VDF3,fontsize=22,weight='bold')
    ax_VDF1_per.text(-1.5,-5.0,'t = '+str(step/2.0)+'s',fontsize=22,weight='bold')
    figname=run+"_VDFs_"+str(step).rjust(7,'0')
    extension = '.png'
    plt.savefig(outputLocation+figname+extension,dpi=300)
    print(outputLocation+figname+extension)
    plt.close()


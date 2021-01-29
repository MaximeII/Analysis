import pytools as pt
import sys, os, socket
import numpy as np
import matplotlib.pyplot as plt
import colormaps as cmaps
import matplotlib
import matplotlib.patches as patches
from multiprocessing import Pool


# Avoids opening a figure window
if str(matplotlib.get_backend()) != 'Agg':
    plt.switch_backend('Agg') 

numproc   = 20
run       = sys.argv[1]
bulkStart = int(sys.argv[2])
bulkEnd   = int(sys.argv[3])

scale_map = 2.5
scale_vdf = 2.5

RE = 6371e3

#fileLocation='/wrk/users/dubart/'+run+'km/maxwellian_periodic/'
fileLocation='/wrk/users/dubart/'+run+'km/diffusion/bulk/'

panel = '$(a)$'
panel_VDF1 = '$(b)$'
panel_VDF2 = '$(c)$'
panel_VDF3 = '$(d)$'

outputLocation=outputdir='/wrk/users/dubart/analysis/subgrid/300_test/diff_test/fig/Colormap/'
path_save = '/wrk/users/dubart/analysis/subgrid/300_test/diff_test/data/'


def plot_cmap(step):

    # Source data file
    bulkname = "bulk."+str(step).rjust(7,'0')+".vlsv"
    print(bulkname)

    # CellIDs of interest
    cid1 = int(np.load(path_save+'CellID_'+run+'.npy'))
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
    pt.plot.plot_colormap(filename=fileLocation+bulkname,
                          var='vg_b_vol',
                          vmin=2.2e-8,vmax=2.7e-8,
                          boxm=[coords[0]-2.0*RE,coords[0]+2.0*RE,coords[1]-2.0*RE,coords[1]+2.0*RE],
                          run=run,
                          colormap='hot_desaturated',
                          scale=scale_map,
                          axes=ax_col,cbaxes=cax_col,
                          outputdir=outputLocation,lin=1)

    # Addition of location of cells with VDFs
    vlsvReader=pt.vlsvfile.VlsvReader(fileLocation+bulkname)
    xCid,yCid,zCid = vlsvReader.get_cell_coordinates(cid1)
    VDFcoord_1 = [xCid/Re,yCid/Re,zCid/Re]

    ax_col.scatter(VDFcoord_1[0], VDFcoord_1[1], color='black',marker='o',s=150)
    ax_col.scatter(VDFcoord_1[0], VDFcoord_1[1], color='white',marker='o',s=70)

    # VDFs of first cell
    pt.plot.plot_vdf(filename=fileLocation+bulkname,
                     cellids=[cid1],
                     bpara1=1,
                     box=[-2.5e6,2.5e6,-2.5e6,2.5e6],
                     fmin=1e-15, fmax=1e-11,
                     slicethick=1,
                     axes=ax_VDF1_par,cbaxes=cax_VDF,
                     cbulk=1,
                     axisunit=6,
                     title='',
                     scale=scale_vdf,
                     cbtitle='',
                     colormap='nipy_spectral')

    pt.plot.plot_vdf(filename=fileLocation+bulkname,
                     cellids=[cid1],
                     bperp=1,
                     box=[-2.5e6,2.5e6,-2.5e6,2.5e6],
                     fmin=1e-15, fmax=1e-11,
                     slicethick=1,
                     axes=ax_VDF1_per,nocb=1,
                     cbulk=1,
                     axisunit=6,
                     title='',
                     scale=scale_vdf,
                     colormap='nipy_spectral')

    pt.plot.plot_vdf(filename=fileLocation+bulkname,
                     cellids=[cid1],
                     bpara=1,
                     box=[-2.5e6,2.5e6,-2.5e6,2.5e6],
                     fmin=1e-15, fmax=1e-11,
                     slicethick=1,
                     axes=ax_VDF1_per2,nocb=1,
                     cbulk=1,
                     axisunit=6,
                     title='',
                     scale=scale_vdf,
                     colormap='nipy_spectral')

    ax_col.text(-1.3,21.1,panel,fontsize=25,weight='bold')
    ax_VDF1_par.text(10.0,2.8,'$f(v)[m^{-6}s^{3}]$' ,fontsize=20,weight='bold')
    ax_VDF1_par.text(-2.3,2.0,panel_VDF1,fontsize=22,weight='bold')
    ax_VDF1_per.text(-2.3,2.0,panel_VDF2,fontsize=22,weight='bold')
    ax_VDF1_per2.text(-2.3,2.0,panel_VDF3,fontsize=22,weight='bold')
    fig.suptitle(run,fontsize=20)
    figname=run+"_VDFs_diff_"+str(step).rjust(7,'0')
    extension = '.png'
    plt.savefig(outputLocation+figname+extension,dpi=300)
    print(outputLocation+figname+extension)


pool = Pool(numproc)
pool.map(plot_cmap,range(bulkStart,bulkEnd+1))
pool.terminate()

import pytools as pt
import sys, os, socket
import numpy as np
import matplotlib.pyplot as plt
import colormaps as cmaps
import matplotlib
import matplotlib.patches as patches
from multiprocessing import Pool

numproc = 20

scale_map = 2.5
scale_vdf = 2.5

path_bulk = '/wrk/group/spacephysics/vlasiator/2D/BCQ/restart/'
path_save = '/wrk/users/dubart/analysis/subgrid/data/BCQ/restart/'
path_fig  = '/wrk/users/dubart/analysis/subgrid/fig/BCQ/restart/'

bulkname = 'restart.0001361.vlsv'

CellIDs = np.load(path_save+'CellIDs_va.npy')

f = pt.vlsvfile.VlsvReader(path_bulk+bulkname)

Re = 6371e3 

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

def plot_cmap(step):

    fig = plt.figure(figsize=(8.*figratio,8.))
    # Colormap and its colorbar
    ax_col = fig.add_axes([mh,mv,l1,h1])
    cax_col = fig.add_axes([mh+l1+ecb,mv,lcb,h1])
    # VDFs and their colorbar
    deltaX = mh+l1+ecb+lcb+2*eh
    ax_VDF1_par  = fig.add_axes([deltaX,mv+h2+ev,l2,h2])
    ax_VDF1_per  = fig.add_axes([deltaX+l2+eh,mv+h2+ev,l2,h2])
    ax_VDF1_per2 = fig.add_axes([deltaX,mv,l2,h2])
    ax_VDF2_per = fig.add_axes([deltaX+l2+eh,mv,l2,h2])
    cax_VDF = fig.add_axes([deltaX+2*l2+eh+1.5*ecb,mv,lcb,h1-0.01])
    
    
    # Colormap with fg_b
    pt.plot.plot_colormap(filename=path_bulk+bulkname,
                          var='b',
                          vmin=1.0e-9,vmax=3.0e-8,
                          boxre=[-40.0,20.0,-20.0,30.0],
                          run='BCQ',
                          colormap='hot_desaturated',
                          scale=scale_map,
                          axes=ax_col,cbaxes=cax_col,
                          outputdir=path_fig,lin=1)
    
    # Addition of location of cells with VDFs
    vlsvReader=pt.vlsvfile.VlsvReader(path_bulk+bulkname)
    xCid,yCid,zCid = vlsvReader.get_cell_coordinates(CellIDs[step])
    VDFcoord_1 = [xCid/Re,yCid/Re,zCid/Re]
    
    ax_col.scatter(VDFcoord_1[0], VDFcoord_1[2], color='black',marker='o',s=150)
    ax_col.scatter(VDFcoord_1[0], VDFcoord_1[2], color='white',marker='o',s=70)
    
    # VDFs of first cell
    pt.plot.plot_vdf(filename=path_bulk+bulkname,
                     cellids=CellIDs[step],
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
    
    pt.plot.plot_vdf(filename=path_bulk+bulkname,
                     cellids=CellIDs[step],
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
    
    pt.plot.plot_vdf(filename=path_bulk+bulkname,
                     cellids=CellIDs[step],
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
    
    ax_VDF1_par.text(10.0,2.8,'$f(v)[m^{-6}s^{3}]$' ,fontsize=20,weight='bold')
    plt.savefig(path_fig+'VDFs_restart_'+str(step).rjust(4,'0')+'.png',dpi=300)
    print(path_fig+'VDFs_restart_'+str(step).rjust(4,'0')+'.png')

    return

pool = Pool(numproc)
pool.map(plot_cmap,range(0,len(CellIDs)))
pool.terminate()


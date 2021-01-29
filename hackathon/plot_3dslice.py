import pytools as pt
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys

run  = sys.argv[1]
bulk = sys.argv[2]

path_fig  = '/wrk/users/dubart/analysis/hackathon/fig/'
path_bulk = '/wrk/group/spacephysics/vlasiator/3D/'+run+'/bulk/dense_cold_hall1e5_afterRestart374/'

filename = 'bulk1.'+str(bulk).rjust(7,'0')+'.vlsv'

pt.plot.plot_colormap3dslice(filename = path_bulk+filename,
                             outputdir = path_fig,
                             vmin = -1.0e-8,vmax = 1.0e-8,
                             #title ='',
                             boxre=[0.0,15.0,-25.0,25.0],
                             var = 'vg_b_vol', operator ='z',
                             run = run)

                             

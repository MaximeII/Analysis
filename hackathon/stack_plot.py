import pytools as pt
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys
import matplotlib.colors as colors
import scipy
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


bulks = np.arange(662,1470,1)        

path_save = '/wrk/users/dubart/analysis/hackathon/data/'
path_fig  = '/wrk/users/dubart/analysis/hackathon/fig/'

components = ['Bmag','Bx','By','Bz','np','Vpx','Vpy','Vpz','Tpara','Tperp','dBzdx']
units      = ['(T)','($m^{-3}$)','(m/s)','(K)','(T/m)']
vmins      = [-3.0e-8,0.0e6,-1.5e6,10e6,-1.0e-15]
vmaxs      = [3.0e-8,2.5e6,1.5e6,100e6,1.0e-15]

filetitle = 'EGI_cut_x_y0_z0'
    
RE = 6371e3

# Read the first file to determine the number of lines and of parameters in the file

input_f = open(path_save+filetitle+'_'+str(bulks[0]).rjust(7,'0')+".dat", 'r')
    
with input_f as f:
    header = f.readline().split()
    lines = (line for line in f if not line.startswith('x'))
    alldata = np.loadtxt(lines, skiprows=1)

nlines = alldata.shape[0]
line_coordinates = np.array((nlines,3))
line_coordinates = alldata[:,0:2]

ntimes = len(bulks)

# Create a dictionary that will contain all the data
data_dict = {}

for param in header:
    data_dict[param] = np.empty((ntimes,nlines))


#read in data
for i in range(len(bulks)):
    input_f = open(path_save+filetitle+'_'+str(bulks[i]).rjust(7,'0')+".dat", 'r')
    with input_f as f:
        lines = (line for line in f if not line.startswith('x'))
        alldata = np.loadtxt(lines, skiprows=1)
        #alldata = np.transpose(alldata)
        
        for param,ind in zip(header,range(0,len(header))):
            data_dict[param][i,:] = alldata[:,ind]

Bz = np.transpose(data_dict['Bz'])
dx = abs(line_coordinates[0,0] - line_coordinates[1,0])

dBzdx = np.zeros([Bz.shape[0],Bz.shape[1]])

dBzdx[0]    = (Bz[1,:] - Bz[0,:]) / dx
dBzdx[1:-1] = (Bz[2:,:] - Bz[0:-2,:]) / (2*dx)
dBzdx[-1]   = (Bz[-1,:] - Bz[-2,:] / dx)      

for i in range(0,len(components)):
 
    if components[i][0] == 'd':
        unit = units[4]
        vmin = vmins[4]
        vmax = vmaxs[4]
        cmapname = 'seismic'
        parameter = dBzdx
    else:
        parameter = np.transpose(data_dict[components[i]])

    if components[i][0] == 'B':
        unit = units[0]
        vmin = vmins[0]
        vmax = vmaxs[0]
        cmapname = 'seismic'

    elif components[i][0] == 'n':
        unit = units[1]
        vmin = vmins[1]
        vmax = vmaxs[1]
        cmapname = 'viridis'


    elif components[i][0] == 'V':
        unit = units[2]
        vmin = vmins[2]
        vmax = vmaxs[2]
        cmapname = 'seismic'

    elif components[i][0] == 'T':
        unit = units[3]
        vmin = vmins[3]
        vmax = vmaxs[3]
        cmapname = 'viridis'

    if components[i] == 'Bmag':
        vmin = 0.0
        cmapname = 'viridis' 


    cmapuse = matplotlib.cm.get_cmap(name=cmapname)
    #norm  = colors.LogNorm(vmin=-3.0e-9, vmax=3.0e-8)

    fig, ax = plt.subplots()

    [XmeshXY,YmeshXY] = np.meshgrid(np.linspace(bulks[0],bulks[-1],num=ntimes),np.linspace(line_coordinates[0,0]/RE,line_coordinates[-1,0]/RE,num=nlines))

    im = ax.pcolormesh(XmeshXY,YmeshXY,parameter, cmap=cmapuse,vmin=vmin,vmax=vmax)
    cbar = fig.colorbar(im,ax=ax)
    if components[i][0] == 'd':
        cbar.set_label('$\\partial B_z / \\partial x$'+' '+unit)
    else:
        cbar.set_label(components[i]+' '+unit)

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    plt.grid(True,which='both')

    plt.xlabel('Time (s)')
    plt.ylabel('X (RE)')
    plt.title('Y = -4.0, Z = 0.0')


    filetitle = 'EGI_cut_x_y-4_z0'
    plt.savefig(path_fig+filetitle+'_'+components[i]+'.png',dpi=300)
    print(path_fig+filetitle+'_'+components[i]+'.png')
    plt.close()

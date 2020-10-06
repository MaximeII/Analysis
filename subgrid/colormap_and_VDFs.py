cpmport pytools as pt
import sys, os, socket
import numpy as np
import matplotlib.pyplot as plt
import colormaps as cmaps
import matplotlib
import matplotlib.patches as patches

# Avoids opening a figure window
if str(matplotlib.get_backend()) is not 'Agg':
    plt.switch_backend('Agg') 

run = sys.argv[1]
scale_map = 2.5
scale_vdf = 2.5

if run == 'BGA':
    fileLocation='/scratch/project_2000203/2D/BGA/reverted_ionosphere_field_boundary/'
    panel = '$(e)$'
    panel_VDF1 = '$(f)$'
    panel_VDF2 = '$(g)$'
    panel_VDF3 = '$(h)$'
    reso  = ' $600$'
    var = 'fg_b'
elif run == 'BCQ':
    fileLocation='/scratch/project_2000203/2D/'+run+'/bulk/'
    panel = '$(a)$'
    panel_VDF1 = '$(b)$'
    panel_VDF2 = '$(c)$'
    panel_VDF3 = '$(d)$'
    reso = ' $300$'
    var = 'b'
elif run == 'BCG':
    fileLocation='/scratch/project_2000203/2D/'+run+'/bulk/'
    panel = '$(i)$'
    panel_VDF1 = '$(j)$'
    panel_VDF2 = '$(k)$'
    panel_VDF3 = '$(l)$'
    reso = ' $900$'
    var = 'b'

outputLocation=outputdir='/users/dubartma/analysator/Data_Analysis/'+run+'/Fig/'
#path_save = '/homeappl/home/dubartma/appl_taito/analysator/Whistler/ABA/Data/'
path_save = '/users/dubartma/analysator/Data_Analysis/'+run+'/Data/'

# Frame extent for this job given as command-line arguments
if len(sys.argv)==4: # Starting and end frames given
    timetot = range(int(sys.argv[2]), int(sys.argv[3]), 1)
else: # Only starting frame given, generate one frame
    timetot = range(int(sys.argv[2]), int(sys.argv[2])+1, 1)
for j in timetot:
    # Source data file
    bulkname = "bulk."+str(j).rjust(7,'0')+".vlsv"
    print(bulkname)

    #f = pt.vlsvfile.VlsvReader(fileLocation+bulkname)
    #
    ##read in data B
    #B = f.read_fsgrid_variable("fg_b")[:,0,:]
    #B_mag = numpy.linalg.norm(B,axis=2).T

    # CellIDs of interest
    cid1 = int(np.load(path_save+'CellID_SC2.npy'))
    Re = 6371e3 # m

    f = pt.vlsvfile.VlsvReader(fileLocation+bulkname)
    coords = f.get_cell_coordinates(cid1)

    if run == 'BGA':
        B      = f.read_fsgrid_variable('fg_b')
        v      = f.read_variable("vg_v",cid1)

        xsize  = f.read_parameter('xcells_ini')
        xmax   = f.read_parameter("xmax")
        xmin   = f.read_parameter("xmin")
        dx     = (xmax-xmin)/xsize

        zmax   = f.read_parameter("zmax")
        zmin   = f.read_parameter("zmin")

        zero_x = int(- xmin/dx)
        zero_z = int(- zmin/dx)

        coord = coords/dx
        B     = B[int(coord[0]-0.5 + zero_x),0,int(coord[2]-0.5 + zero_z)]   
    else:
        B = f.read_variable('B',cid1)
        v = f.read_variable("v",cid1)

    Bnorm = np.linalg.norm(B)
    vnorm = np.linalg.norm(v)

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
                          var=var,
                          vmin=1.0e-8,vmax=3.0e-8,
                          boxre=[0.0,9.0,12.0,21.0],
                          run=run,
                          colormap='hot_desaturated',
                          scale=scale_map,
                          axes=ax_col,cbaxes=cax_col,
                          outputdir=outputLocation,lin=1)

    ax_col.arrow(coords[0]/Re,coords[2]/Re-0.5,B[0]/Bnorm,B[2]/Bnorm,color='white',width=0.1)
    ax_col.arrow(coords[0]/Re,coords[2]/Re+0.5,v[0]/vnorm,v[2]/vnorm,color='black',width=0.1)
    ax_col.text(coords[0]/Re-0.5,coords[2]/Re-1.5,'$\\vec{B}$',color='white',fontsize=35,weight='bold')
    ax_col.text(coords[0]/Re-1.2,coords[2]/Re+1.5,'$\\vec{v}$',color='black',fontsize=35,weight='bold')

    #ax_col.pcolormesh(B_mag,vmin=1.0e-8,vmax=3.0e-8,cmap='hot_desaturated')
    #ax_col.set_xlabel('$X[R_E]$')
    #ax_col.set_ylabel('$Z[R_E]$')
    #ax_col.set_title('$t = '+str(j/2)+' s')    

    rect = patches.Rectangle((3.0,15.0),3.0,3.0,linewidth=2,edgecolor='red',facecolor='none')
    ax_col.add_patch(rect)

    # Addition of location of cells with VDFs
    vlsvReader=pt.vlsvfile.VlsvReader(fileLocation+bulkname)
    xCid,yCid,zCid = vlsvReader.get_cell_coordinates(cid1)
    VDFcoord_1 = [xCid/Re,yCid/Re,zCid/Re]

    ax_col.scatter(VDFcoord_1[0], VDFcoord_1[2], color='black',marker='o',s=150)
    ax_col.scatter(VDFcoord_1[0], VDFcoord_1[2], color='white',marker='o',s=70)

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

   # # VDFs of second cell
   # pt.plot.plot_vdf(filename=fileLocation+bulkname,
   #                  cellids=[cid2],
   #                  bpara1=1,
   #                  box=[-2.5e6,2.5e6,-2.5e6,2.5e6],
   #                  fmin=1e-15, fmax=1e-11,
   #                  slicethick=1,
   #                  axes=ax_VDF2_par,nocb=1,
   #                  cbulk=1,
   #                  axisunit=6,
   #                  title='',
   #                  scale=1.15,
   #                  colormap='nipy_spectral')

   # pt.plot.plot_vdf(filename=fileLocation+bulkname,
   #                  cellids=[cid2],
   #                  bperp=1,
   #                  box=[-2.5e6,2.5e6,-2.5e6,2.5e6],
   #                  fmin=1e-15, fmax=1e-11,
   #                  slicethick=1,
   #                  axes=ax_VDF2_per,nocb=1,
   #                  cbulk=1,
   #                  axisunit=6,
   #                  title='',
   #                  scale=1.15,
   #                  colormap='nipy_spectral')
    ax_col.text(-1.3,21.1,panel,fontsize=25,weight='bold')
    ax_VDF1_par.text(10.0,2.8,'$f(v)[m^{-6}s^{3}]$' ,fontsize=20,weight='bold')
    ax_VDF1_par.text(-2.3,2.0,panel_VDF1,fontsize=22,weight='bold')
    ax_VDF1_per.text(-2.3,2.0,panel_VDF2,fontsize=22,weight='bold')
    ax_VDF1_per2.text(-2.3,2.0,panel_VDF3,fontsize=22,weight='bold')
    fig.suptitle('$\\Delta r = $'+reso+' $km$',fontsize=20)
    figname=run+"_VDFs_"+str(j).rjust(7,'0')
    extension = '.png'
    plt.savefig(outputLocation+figname+extension,dpi=300)
    print(outputLocation+figname+extension)
    #plt.savefig(outputLocation+figname+'.eps')

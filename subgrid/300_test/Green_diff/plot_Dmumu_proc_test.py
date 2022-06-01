import pytools as pt
import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt

run    = 'proc_test'

path_fig  = '/wrk-vakka/users/dubart/analysis/subgrid/300_test/Green_diff/fig/mu/'+run+'/'
path_save = '/wrk/users/dubart/analysis/subgrid/300_test/Green_diff/data/mu/'+run+'/'
path_bulk = '/wrk/users/dubart/diff_test/'+run+'/mu_bulk/'

# Figure layout
xmargin   = 0.10
ymargin   = 0.11
vdfwidth  = 0.25
vdfheight = 0.35
axwidth   = 0.45
axheight  = 0.39
cbmargin  = 0.015
cbwidth   = 0.02
axspace   = 0.06

if len(sys.argv) == 3: # Starting and end frames given
    timetot = range(int(sys.argv[1]), int(sys.argv[2]),1)
else: # Only starting frame given, generate one frame
    timetot = range(int(sys.argv[1]), int(sys.argv[1])+1, 1)
for step in timetot:

    fig = plt.figure(figsize=(15,8))

    if step == 0:
        bulkname       = "initial-grid.0000000.vlsv"
    else:
        bulkname = "bulk."+str(step).rjust(7,'0')+".vlsv"

    ax_fmu      = fig.add_axes([xmargin,ymargin+axheight+axspace*2-0.08,axwidth,axheight])
    ax_dmumu    = fig.add_axes([xmargin,ymargin,axwidth,axheight])
    ax_vdf_300  = fig.add_axes([xmargin+axwidth+axspace*1.5,ymargin+vdfheight+axspace*2,vdfwidth,vdfheight])
    #ax_vdf_diff = fig.add_axes([xmargin+axwidth+axspace*1.5,ymargin,vdfwidth,vdfheight])
    vdf_cax     = fig.add_axes([xmargin+axwidth+axspace*1.5+vdfwidth+cbmargin,ymargin,cbwidth,vdfheight*2+axspace*2])

    #ax_vdf_diff.xaxis.label.set_color('red')
    #ax_vdf_diff.yaxis.label.set_color('red')
    #ax_vdf_diff.tick_params(axis='both',colors='red')
    #ax_vdf_diff.spines['top'].set_color('red')
    #ax_vdf_diff.spines['bottom'].set_color('red')
    #ax_vdf_diff.spines['left'].set_color('red')
    #ax_vdf_diff.spines['right'].set_color('red')

    ax_fmu.set_xlim(-1,1)
    ax_fmu.set_ylim(0.1,1e6)
    ax_fmu.set_yscale('log')
    ax_fmu.set_ylabel('n [m$^{-3}$]',fontsize=20,weight='bold')    
    
    ax_fmu.xaxis.set_ticklabels([])
    ax_fmu.xaxis.set_tick_params(labelsize=18)
    ax_fmu.yaxis.set_tick_params(labelsize=19)

    ax_dmumu.set_xlim(-1,1)
    ax_dmumu.set_ylim(0.0,0.3)
    ax_dmumu.set_xlabel('$\\mu$',fontsize=20,weight='bold')
    ax_dmumu.set_ylabel('|D$_{\\mu\\mu}$| [s$^{-1}$]',fontsize=20,weight='bold')

    ax_dmumu.xaxis.set_tick_params(labelsize=17)
    ax_dmumu.yaxis.set_tick_params(labelsize=19)

    labels = ax_fmu.get_xticklabels() + ax_fmu.get_yticklabels() + ax_dmumu.get_xticklabels() + ax_dmumu.get_yticklabels()
    [label.set_fontweight('bold') for label in labels]

    CellID = 1

    # plot fmu and Dmumu
    file_mu = open(path_save + 'mu_'+run+'_'+str(step).rjust(7,'0')+'.txt','r')
    Lines   = file_mu.readlines()
    file_mu.close()

    mu  = np.fromstring(Lines[3],dtype='float',sep=' ')
    fmu = np.fromstring(Lines[5],dtype='float',sep=' ')

    mu_mid = (mu[:-1] + mu[1:])/2

    ax_fmu.plot(mu_mid,fmu,linewidth=2,color='black')

    Dmumu = np.load(path_save+'Dmumu_'+run+'.npy')

    ax_dmumu.plot(mu_mid,abs(Dmumu[:,step]),linewidth=2,color='black')

    ax_dmumu.axhline(y=0.01,color='red',linestyle='--')

    #try:
    # plot VDFs

    pt.plot.plot_vdf(filename   = path_bulk+bulkname,
                     cellids    = [CellID],
                     xy         = 1,
                     box        = [-2.0e6,2.0e6,-2.0e6,2.0e6],
                     fmin       = 1e-15, fmax = 1e-11,
                     cbulk      = 1,
                     axes       = ax_vdf_300,cbaxes=vdf_cax,
                     slicethick = 1,
                     axisunit   = 6,
                     title      = '',
                     scale      = 2.5,
                     cbtitle    = '',
                     #setThreshold = 1e-21,
                     colormap   = 'viridis')
    
    #except:
    #    print("Something went wrong")


    cb_text = r'$f($v$)$ [m$^{-6}$ s$^{3}$]'
    ax_vdf_300.text(-1.9,1.6,'(c)',fontsize = 22, weight ='bold')

    ax_fmu.text(0.01,0.9,'(a)',transform=ax_fmu.transAxes,fontsize=22,weight='bold')
    ax_dmumu.text(0.01,0.9,'(b)',transform=ax_dmumu.transAxes,fontsize=22,weight='bold')

    plt.suptitle('t = '+str(step/100)+'s',fontsize=18,weight='bold')
    plt.savefig(path_fig+'Dmumu_'+run+'_'+str(step).rjust(7,'0')+'.png',dpi=300)
    print(path_fig+'Dmumu_'+run+'_'+str(step).rjust(7,'0')+'.png')
    plt.close()

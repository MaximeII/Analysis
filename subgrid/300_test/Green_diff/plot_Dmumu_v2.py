import pytools as pt
import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt

run = sys.argv[1]

runs    = ['300',run]
colours = ['black','red']
legends = ['300','900+diff']

path_fig  = '/wrk-vakka/users/dubart/analysis/subgrid/300_test/Green_diff/fig/mu/'

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

bulkStart = int(sys.argv[2])
bulkEnd   = int(sys.argv[3])
timetot   = bulkEnd - bulkStart +1

tsize = 4 #timetot

for step in range(0,tsize):

    fig = plt.figure(figsize=(15,8))

    if step == 0:
        bulkname = "initial-grid.0000000.vlsv"
    else:
        bulkname = "bulk."+str(step*100).rjust(7,'0')+".vlsv"

    ax_fmu      = fig.add_axes([xmargin,ymargin+axheight+axspace*2-0.08,axwidth,axheight])
    ax_dmumu    = fig.add_axes([xmargin,ymargin,axwidth,axheight])
    ax_vdf_300  = fig.add_axes([xmargin+axwidth+axspace*1.5,ymargin+vdfheight+axspace*2,vdfwidth,vdfheight])
    ax_vdf_diff = fig.add_axes([xmargin+axwidth+axspace*1.5,ymargin,vdfwidth,vdfheight])
    vdf_cax     = fig.add_axes([xmargin+axwidth+axspace*1.5+vdfwidth+cbmargin,ymargin,cbwidth,vdfheight*2+axspace*2])

    ax_vdf_diff.xaxis.label.set_color(colours[1])
    ax_vdf_diff.yaxis.label.set_color(colours[1])
    ax_vdf_diff.tick_params(axis='both',colors=colours[1])
    ax_vdf_diff.spines['top'].set_color(colours[1])
    ax_vdf_diff.spines['bottom'].set_color(colours[1])
    ax_vdf_diff.spines['left'].set_color(colours[1])
    ax_vdf_diff.spines['right'].set_color(colours[1])

    ax_fmu.set_xlim(-1,1)
    #ax_fmu.set_ylim(0.1,1e6)
    ax_fmu.set_ylim(0.1,1e6)
    #ax_fmu.set_yscale('log')
    ax_fmu.set_ylabel('n [m$^{-3}$]',fontsize=20,weight='bold')    
    
    ax_fmu.xaxis.set_ticklabels([])
    ax_fmu.xaxis.set_tick_params(labelsize=18)
    ax_fmu.yaxis.set_tick_params(labelsize=19)

    ax_dmumu.set_xlim(-1,1)
    ax_dmumu.set_ylim(-0.01,0.01)
    ax_dmumu.set_xlabel('$\\mu$',fontsize=20,weight='bold')
    ax_dmumu.set_ylabel('D$_{\\mu\\mu}$ [s$^{-1}$]',fontsize=20,weight='bold')

    ax_dmumu.xaxis.set_tick_params(labelsize=17)
    ax_dmumu.yaxis.set_tick_params(labelsize=19)

    labels = ax_fmu.get_xticklabels() + ax_fmu.get_yticklabels() + ax_dmumu.get_xticklabels() + ax_dmumu.get_yticklabels()
    [label.set_fontweight('bold') for label in labels]

    for i in range(0,len(runs)):

        if runs[i] == '300':
            path_bulk = '/wrk-vakka/users/dubart/diff_test/dmumu/'+runs[i]+'/bulk/'
        else:
            path_bulk = '/wrk-vakka/users/dubart/diff_test/dmumu/diff/'+runs[i]+'/CFL1.0/bulk/'
            #path_bulk = '/wrk-vakka/users/dubart/diff_test/dmumu/diff/'+runs[i]+'/opti_v1/bulk/'

        path_save = '/wrk-vakka/users/dubart/analysis/subgrid/300_test/Green_diff/data/mu/'+runs[i]+'/'
        CellID = int(np.load(path_save+'CellID_'+runs[i]+'.npy'))

        # plot fmu and Dmumu
        file_mu = open(path_save + 'mu_'+runs[i]+'_'+str(step*100).rjust(7,'0')+'.txt','r')
        Lines   = file_mu.readlines()
        file_mu.close()

        mu  = np.fromstring(Lines[3],dtype='float',sep=' ')
        fmu = np.log10(np.fromstring(Lines[5],dtype='float',sep=' '))
 
        mu_mid = (mu[:-1] + mu[1:])/2

        ax_fmu.plot(mu_mid,fmu,linewidth=2,color=colours[i],label=legends[i])

        Dmumu = np.load(path_save+'Dmumu_'+runs[i]+'.npy')

        ax_dmumu.plot(mu_mid,Dmumu[:,step],linewidth=2,color=colours[i])
        #ax_dmumu.scatter(mu_mid,Dmumu[:,step],color=colours[i],s=10)

        ax_dmumu.axhline(y=0.01,color='red',linestyle='--')

        #try:
        # plot VDFs
        if i == 0:

            pt.plot.plot_vdf(filename   = path_bulk+bulkname,
                             cellids    = [CellID],
                             xy         = 1,
                             box        = [-2.0e6,2.0e6,-2.0e6,2.0e6],
                             fmin       = 1e-15, fmax = 1e-11,
                             cbulk      = 1,
                             axes       = ax_vdf_300,nocb=1,
                             slicethick = 1,
                             axisunit   = 6,
                             title      = '',
                             scale      = 2.5,
                             cbtitle    = '',
                             #setThreshold = 1e-21,
                             colormap   = 'viridis')
       
        elif i == 1:
       
            pt.plot.plot_vdf(filename   = path_bulk+bulkname,
                             cellids    = [CellID],
                             xy         = 1,
                             box        = [-2.0e6,2.0e6,-2.0e6,2.0e6],
                             fmin       = 1e-15, fmax = 1e-11,
                             cbulk      = 1,
                             axes       = ax_vdf_diff,cbaxes=vdf_cax,
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
    ax_vdf_diff.text(2.0,-3.0,cb_text,fontsize=20,weight='bold')

    ax_fmu.legend(loc = 'lower center',fontsize=20)

    ax_vdf_300.text(-1.9,1.6,'(c)',fontsize = 22, weight ='bold')
    ax_vdf_diff.text(-1.9,1.6,'(d)',fontsize = 22, weight ='bold')

    ax_fmu.text(0.01,0.9,'(a)',transform=ax_fmu.transAxes,fontsize=22,weight='bold')
    ax_dmumu.text(0.01,0.9,'(b)',transform=ax_dmumu.transAxes,fontsize=22,weight='bold')

    plt.suptitle('t = '+str(round(step*100,1))+'s',fontsize=18,weight='bold')
    plt.savefig(path_fig+'Dmumu_'+str(step).rjust(7,'0')+'.png',dpi=300)
    print(path_fig+'Dmumu_'+str(step).rjust(7,'0')+'.png')
    plt.close()

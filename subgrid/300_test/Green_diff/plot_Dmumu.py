import pytools as pt
import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt

run = sys.argv[1]

runs    = ['300',run]
colours = ['black','red']
legends = ['300','900+diff']

bulkStart = sys.argv[2]
bulkEnd   = sys.argv[3]

path_fig  = '/wrk-vakka/users/dubart/analysis/subgrid/300_test/Green_diff/fig/mu/'

dt = 0.1

Deltat = (int(bulkEnd) - int(bulkStart))*0.1

# Figure layout
xmargin   = 0.08
ymargin   = 0.11
vdfwidth  = 0.25
vdfheight = 0.35
axwidth   = 0.33
axheight  = 0.39
cbmargin  = 0.015
cbwidth   = 0.02
axspacew  = 0.04
axspaceh  = 0.06

fig = plt.figure(figsize=(20,8))

if bulkStart == '0':
    bulkname_start = "initial-grid.0000000.vlsv"
else:
    bulkname_start = "bulk."+str(bulkStart).rjust(7,'0')+".vlsv"
bulkname_end = "bulk."+str(bulkEnd).rjust(7,'0')+".vlsv"

ax_fmu            = fig.add_axes([xmargin,ymargin+axheight+axspaceh*2-0.08,axwidth,axheight])
ax_dmumu          = fig.add_axes([xmargin,ymargin,axwidth,axheight])
ax_vdf_300_start  = fig.add_axes([xmargin+axwidth+axspacew,ymargin+vdfheight+axspaceh*2,vdfwidth,vdfheight])
ax_vdf_300_end    = fig.add_axes([xmargin+axwidth+axspacew*0.8+vdfwidth,ymargin+vdfheight+axspaceh*2,vdfwidth,vdfheight])
ax_vdf_diff_start = fig.add_axes([xmargin+axwidth+axspacew,ymargin,vdfwidth,vdfheight])
ax_vdf_diff_end   = fig.add_axes([xmargin+axwidth+axspacew*0.8+vdfwidth,ymargin,vdfwidth,vdfheight])
vdf_cax     = fig.add_axes([xmargin+axwidth+axspacew*0.3+vdfwidth*2,ymargin,cbwidth,vdfheight*2+axspaceh*2])

ax_vdf_diff_start.xaxis.label.set_color(colours[1])
ax_vdf_diff_start.yaxis.label.set_color(colours[1])
ax_vdf_diff_start.tick_params(axis='both',colors=colours[1])
ax_vdf_diff_start.spines['top'].set_color(colours[1])
ax_vdf_diff_start.spines['bottom'].set_color(colours[1])
ax_vdf_diff_start.spines['left'].set_color(colours[1])
ax_vdf_diff_start.spines['right'].set_color(colours[1])

ax_vdf_diff_end.xaxis.label.set_color(colours[1])
ax_vdf_diff_end.yaxis.label.set_color(colours[1])
ax_vdf_diff_end.tick_params(axis='both',colors=colours[1])
ax_vdf_diff_end.spines['top'].set_color(colours[1])
ax_vdf_diff_end.spines['bottom'].set_color(colours[1])
ax_vdf_diff_end.spines['left'].set_color(colours[1])
ax_vdf_diff_end.spines['right'].set_color(colours[1])

ax_fmu.set_xlim(-1,1)
ax_fmu.set_ylim(0.1,2.5e5)
#ax_fmu.set_ylim(3,6)
#ax_fmu.set_yscale('log')
ax_fmu.set_ylabel('n [m$^{-3}$]',fontsize=20,weight='bold')    
ax_fmu.ticklabel_format(axis='y',style='sci',scilimits=(5,5))

ax_fmu.xaxis.set_ticklabels([])
ax_fmu.xaxis.set_tick_params(labelsize=18)
ax_fmu.yaxis.set_tick_params(labelsize=19)

ax_dmumu.set_xlim(-1,1)
ax_dmumu.set_ylim(0.0,0.01)
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
    file_mu = open(path_save + 'mu_'+runs[i]+'_'+str(bulkStart).rjust(7,'0')+'.txt','r')
    Lines   = file_mu.readlines()
    file_mu.close()

    mu        = np.fromstring(Lines[3],dtype='float',sep=' ')
    fmu_start = np.fromstring(Lines[5],dtype='float',sep=' ')

    mu_mid = (mu[:-1] + mu[1:])/2

    file_mu = open(path_save + 'mu_'+runs[i]+'_'+str(bulkEnd).rjust(7,'0')+'.txt','r')
    Lines   = file_mu.readlines()
    file_mu.close()

    mu        = np.fromstring(Lines[3],dtype='float',sep=' ')
    fmu_end = np.fromstring(Lines[5],dtype='float',sep=' ')

    mu_mid = (mu[:-1] + mu[1:])/2

    fmu = (fmu_start + fmu_end)/2.0

    ax_fmu.plot(mu_mid,fmu,linewidth=2,color=colours[i],label=legends[i])
 
    Dmumu = np.load(path_save+'Dmumu_integration_'+runs[i]+'_'+str(bulkStart)+'-'+str(bulkEnd)+'.npy')

    ax_dmumu.plot(mu_mid,Dmumu,linewidth=2,color=colours[i])
    #ax_dmumu.scatter(mu_mid,Dmumu[:,step],color=colours[i],s=10)

    ax_dmumu.axhline(y=0.01,color='red',linestyle='--')

    #try:
    # plot VDFs
    if i == 0:

        pt.plot.plot_vdf(filename   = path_bulk+bulkname_start,
                         cellids    = [CellID],
                         xy         = 1,
                         box        = [-2.0e6,2.0e6,-2.0e6,2.0e6],
                         fmin       = 1e-15, fmax = 1e-11,
                         cbulk      = 1,
                         axes       = ax_vdf_300_start,nocb=1,
                         slicethick = 1,
                         axisunit   = 6,
                         title      = '',
                         scale      = 2.5,
                         cbtitle    = '',
                         #setThreshold = 1e-21,
                         colormap   = 'viridis')

        pt.plot.plot_vdf(filename   = path_bulk+bulkname_end,
                         cellids    = [CellID],
                         xy         = 1,
                         box        = [-2.0e6,2.0e6,-2.0e6,2.0e6],
                         fmin       = 1e-15, fmax = 1e-11,
                         cbulk      = 1,
                         axes       = ax_vdf_300_end,nocb=1,
                         slicethick = 1,
                         axisunit   = 6,
                         title      = '',
                         scale      = 2.5,
                         cbtitle    = '',
                         #setThreshold = 1e-21,
                         colormap   = 'viridis')
   
    elif i == 1:
   
        pt.plot.plot_vdf(filename   = path_bulk+bulkname_start,
                         cellids    = [CellID],
                         xy         = 1,
                         box        = [-2.0e6,2.0e6,-2.0e6,2.0e6],
                         fmin       = 1e-15, fmax = 1e-11,
                         cbulk      = 1,
                         axes       = ax_vdf_diff_start,nocb=1,
                         slicethick = 1,
                         axisunit   = 6,
                         title      = '',
                         scale      = 2.5,
                         cbtitle    = '',
                         #setThreshold = 1e-21,
                         colormap   = 'viridis')

        pt.plot.plot_vdf(filename   = path_bulk+bulkname_end,
                         cellids    = [CellID],
                         xy         = 1,
                         box        = [-2.0e6,2.0e6,-2.0e6,2.0e6],
                         fmin       = 1e-15, fmax = 1e-11,
                         cbulk      = 1,
                         axes       = ax_vdf_diff_end,cbaxes=vdf_cax,
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
ax_vdf_diff_end.text(2.0,-3.0,cb_text,fontsize=20,weight='bold')

ax_fmu.legend(loc = 'lower center',fontsize=20)

ax_vdf_300_start.text(-1.9,1.6,'(c)',fontsize = 22, weight ='bold')
ax_vdf_300_end.text(-1.9,1.6,'(d)',fontsize = 22, weight ='bold')
ax_vdf_diff_start.text(-1.9,1.6,'(e)',fontsize = 22, weight ='bold')
ax_vdf_diff_end.text(-1.9,1.6,'(f)',fontsize = 22, weight ='bold')

ax_fmu.text(0.01,0.9,'(a)',transform=ax_fmu.transAxes,fontsize=22,weight='bold')
ax_dmumu.text(0.01,0.9,'(b)',transform=ax_dmumu.transAxes,fontsize=22,weight='bold')

plt.suptitle('$\\Delta$ t = '+str(round(Deltat,1))+'s',fontsize=18,weight='bold')
plt.savefig(path_fig+'Dmumu_'+str(bulkStart)+'-'+str(bulkEnd)+'.png',dpi=300)
print(path_fig+'Dmumu_'+str(bulkStart)+'-'+str(bulkEnd)+'.png')
plt.close()

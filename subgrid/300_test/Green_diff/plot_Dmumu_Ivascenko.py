import pytools as pt
import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt

run = sys.argv[1]

colours = ['black','red']

bulkStart = int(sys.argv[2])
bulkEnd   = int(sys.argv[3])

interval = int(sys.argv[4])

path_fig  = '/wrk-vakka/users/dubart/analysis/subgrid/300_test/Green_diff/fig/mu/'

# Figure layout
xmargin   = 0.12
ymargin   = 0.11
vdfwidth  = 0.25
vdfheight = 0.35
axwidth   = 0.45
axheight  = 0.39
cbmargin  = 0.015
cbwidth   = 0.02
axspacew  = 0.05
axspaceh  = 0.06

if run == '300':
    path_bulk = '/wrk-vakka/users/dubart/diff_test/dmumu/'+run+'/bulk/'
elif run == 'Markus':
    path_bulk = '/wrk-vakka/users/markusb/weirddiffusion/Losscone_may31/'
elif run == 'Markus_artificial':
    path_bulk = '/wrk-vakka/users/markusb/weirddiffusion/Losscone_may31_artificial/'
else:
    path_bulk = '/wrk-vakka/users/dubart/diff_test/dmumu/diff/'+run+'/CFL1.0/bulk/'
    #path_bulk = '/wrk-vakka/users/dubart/diff_test/dmumu/diff/'+runs[i]+'/opti_v1/bulk/'

timetot = range(int(sys.argv[2]), int(sys.argv[3]),1)

for j in timetot:

    bulks = [j, j + interval] 

    fig = plt.figure(figsize=(15,8))
    
    ax_fmu      = fig.add_axes([xmargin,ymargin+axheight+axspaceh*2-0.08,axwidth,axheight])
    ax_dmumu    = fig.add_axes([xmargin,ymargin,axwidth,axheight])
    ax_vdf_300  = fig.add_axes([xmargin+axwidth+axspacew*1.5,ymargin+vdfheight+axspaceh*2,vdfwidth,vdfheight])
    ax_vdf_diff = fig.add_axes([xmargin+axwidth+axspacew*1.5,ymargin,vdfwidth,vdfheight])
    vdf_cax     = fig.add_axes([xmargin+axwidth+axspacew*1.5+vdfwidth+cbmargin,ymargin,cbwidth,vdfheight*2+axspaceh*2])
    
    ax_fmu.set_xlim(-1,1)
    #ax_fmu.set_ylim(0.1,2.5e5)
    ax_fmu.set_ylim(0.1,80)
    #ax_fmu.set_ylim(3,6)
    #ax_fmu.set_yscale('log')
    ax_fmu.set_ylabel('n [m$^{-3}$]',fontsize=20,weight='bold')    
    
    ax_fmu.xaxis.set_ticklabels([])
    ax_fmu.xaxis.set_tick_params(labelsize=18)
    ax_fmu.yaxis.set_tick_params(labelsize=19)
    
    if run != 'Markus':
        ax_fmu.ticklabel_format(axis='y',style='sci',scilimits=(5,5))
    
    ax_vdf_diff.xaxis.label.set_color(colours[1])
    ax_vdf_diff.yaxis.label.set_color(colours[1])
    ax_vdf_diff.tick_params(axis='both',colors=colours[1])
    ax_vdf_diff.spines['top'].set_color(colours[1])
    ax_vdf_diff.spines['bottom'].set_color(colours[1])
    ax_vdf_diff.spines['left'].set_color(colours[1])
    ax_vdf_diff.spines['right'].set_color(colours[1])
    
    
    ax_dmumu.set_xlim(-1,1)
    ax_dmumu.set_ylim(0.0,0.02)
    ax_dmumu.set_xlabel('$\\mu$',fontsize=20,weight='bold')
    ax_dmumu.set_ylabel('D$_{\\mu\\mu}$ [s$^{-1}$]',fontsize=20,weight='bold')
    
    ax_dmumu.xaxis.set_tick_params(labelsize=17)
    ax_dmumu.yaxis.set_tick_params(labelsize=19)
    
    labels = ax_fmu.get_xticklabels() + ax_fmu.get_yticklabels() + ax_dmumu.get_xticklabels() + ax_dmumu.get_yticklabels()
    [label.set_fontweight('bold') for label in labels]


    if bulks[0] == 0:
        bulkname_start = 'initial-grid.0000000.vlsv'
    else:
        bulkname_start = 'bulk.'+str(bulks[0]).rjust(7,'0')+'.vlsv'
    bulkname_end = 'bulk.'+str(bulks[1]).rjust(7,'0')+'.vlsv'
    
    f_start = pt.vlsvfile.VlsvReader(path_bulk + bulkname_start)
    f_end   = pt.vlsvfile.VlsvReader(path_bulk + bulkname_end)
    
    time_start = f_start.read_parameter('time')
    time_end   = f_end.read_parameter('time')
    
    dt = time_end - time_start
    
    path_save = '/wrk-vakka/users/dubart/analysis/subgrid/300_test/Green_diff/data/mu/'+run+'/'
    
    if run == 'Markus_artificial':
        CellID = 1
    else:
        CellID = int(np.load(path_save+'CellID_'+run+'.npy'))
    
    # plot fmu and Dmumu
    file_mu = open(path_save + 'mu_'+run+'_'+str(bulks[0]).rjust(7,'0')+'.txt','r')
    Lines   = file_mu.readlines()
    file_mu.close()
    
    mu        = np.fromstring(Lines[3],dtype='float',sep=' ')
    fmu_start = np.fromstring(Lines[5],dtype='float',sep=' ')
    
    mu_mid = (mu[:-1] + mu[1:])/2
    
    ax_fmu.plot(mu_mid,fmu_start,linewidth=2,color=colours[0],label='t = '+str(round(time_start,1))+'s')
    
    file_mu = open(path_save + 'mu_'+run+'_'+str(bulks[1]).rjust(7,'0')+'.txt','r')
    Lines   = file_mu.readlines()
    file_mu.close()
    
    mu      = np.fromstring(Lines[3],dtype='float',sep=' ')
    fmu_end = np.fromstring(Lines[5],dtype='float',sep=' ')
    
    mu_mid = (mu[:-1] + mu[1:])/2
    
    ax_fmu.plot(mu_mid,fmu_end,linewidth=2,color=colours[1],label='t = '+str(round(time_end,1))+'s')
    
    fmu_avg = (fmu_start + fmu_end)/2.0
    
    ax_fmu.plot(mu_mid,fmu_avg,linewidth=2,color='green')
    
    Dmumu_matrix      = np.load(path_save+'Dmumu_matrix_'+run+'_'+str(j).rjust(7,'0')+'.npy')
    Dmumu_integration = np.load(path_save+'Dmumu_integration_'+run+'_'+str(j).rjust(7,'0')+'.npy')
    
    #ax_dmumu.plot(mu_mid,Dmumu_matrix,linewidth=2,color='blue',label='Matrix')
    ax_dmumu.plot(mu_mid,Dmumu_integration,linewidth=2,color='orange',label='Integration')
    #ax_dmumu.scatter(mu_mid,Dmumu[:,step],color=colours[i],s=10)
    
    #ax_dmumu.axhline(y=0.01,color='black',linestyle='--')
    
    #try:
    # plot VDFs
    
    if run == 'Markus' or run == 'Markus_artificial':
    
        pt.plot.plot_vdf(filename   = path_bulk+bulkname_start,
                         cellids    = [CellID],
                         xy         = 1,
                         box        = [-2.0e6,2.0e6,-2.0e6,2.0e6],
                         #fmin       = 1e-15, fmax = 1e-11,
                         fmin        = 1e-20, fmax = 1e-15,
                         cbulk      = None,
                         axes       = ax_vdf_300,nocb=1,
                         slicethick = 1,
                         axisunit   = 6,
                         title      = '',
                         scale      = 2.5,
                         cbtitle    = '',
                         pop        ='protonmono',
                         setThreshold = 1e-21,
                         colormap   = 'viridis')
        
        
        pt.plot.plot_vdf(filename   = path_bulk+bulkname_end,
                         cellids    = [CellID],
                         xy         = 1,
                         box        = [-2.0e6,2.0e6,-2.0e6,2.0e6],
                         #fmin       = 1e-15, fmax = 1e-11,
                         fmin        = 1e-20, fmax = 1e-15,
                         cbulk      = None,
                         axes       = ax_vdf_diff,cbaxes=vdf_cax,
                         slicethick = 1,
                         axisunit   = 6,
                         title      = '',
                         scale      = 2.5,
                         cbtitle    = '',
                         pop        ='protonmono',
                         setThreshold = 1e-21,
                         colormap   = 'viridis')
    
    else:
    
        pt.plot.plot_vdf(filename   = path_bulk+bulkname_start,
                         cellids    = [CellID],
                         xy         = 1,
                         box        = [-2.0e6,2.0e6,-2.0e6,2.0e6],
                         fmin       = 1e-15, fmax = 1e-11,
                         #fmin        = 1e-20, fmax = 1e-15,
                         cbulk      = 1,
                         axes       = ax_vdf_300,nocb=1,
                         slicethick = 1,
                         axisunit   = 6,
                         title      = '',
                         scale      = 2.5,
                         cbtitle    = '',
                         #pop        ='protonmono',
                         #setThreshold = 1e-21,
                         colormap   = 'viridis')
        
        
        pt.plot.plot_vdf(filename   = path_bulk+bulkname_end,
                         cellids    = [CellID],
                         xy         = 1,
                         box        = [-2.0e6,2.0e6,-2.0e6,2.0e6],
                         fmin       = 1e-15, fmax = 1e-11,
                         #fmin        = 1e-20, fmax = 1e-15,
                         cbulk      = 1,
                         axes       = ax_vdf_diff,cbaxes=vdf_cax,
                         slicethick = 1,
                         axisunit   = 6,
                         title      = '',
                         scale      = 2.5,
                         cbtitle    = '',
                         #pop        ='protonmono',
                         #setThreshold = 1e-21,
                         colormap   = 'viridis')
    #except:
    #    print("Something went wrong")
    
    cb_text = r'$f($v$)$ [m$^{-6}$ s$^{3}$]'
    ax_vdf_diff.text(2.0,-3.0,cb_text,fontsize=20,weight='bold')
    
    ax_fmu.legend(loc = 'upper center',fontsize=20)
    ax_dmumu.legend(loc = 'lower center',fontsize=20)
    
    ax_vdf_300.text(-1.9,1.6,'(c)',fontsize = 22, weight ='bold')
    ax_vdf_diff.text(-1.9,1.6,'(d)',fontsize = 22, weight ='bold')
    
    ax_fmu.text(0.01,0.9,'(a)',transform=ax_fmu.transAxes,fontsize=22,weight='bold')
    ax_dmumu.text(0.01,0.9,'(b)',transform=ax_dmumu.transAxes,fontsize=22,weight='bold')
    
    plt.suptitle('$\\Delta$ t = '+str(round(dt,1))+'s',fontsize=18,weight='bold')
    plt.savefig(path_fig+'Dmumu_'+run+'_'+str(j).rjust(7,'0')+'.png',dpi=300)
    print(path_fig+'Dmumu_'+run+'_'+str(j).rjust(7,'0')+'.png')
    plt.close()

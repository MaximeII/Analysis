import pytools as pt
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os, sys
from multiprocessing import Pool
import scipy
import scipy.ndimage
import matplotlib.colors as colors
import imageio
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


def get_fields(step):

    ''' Auxiliary function for computing magnetic, electric fields and density in parallel

    :kword step: Bulk number of the current file
    returns:     list of fields

    '''

    f = pt.vlsvfile.VlsvReader(filelist[step-start_global])

    # Read cellids. Cellids are randomly sorted because of multiprocessing
    cellids = f.read_variable('CellID')

    # Read field
    B = f.read_variable('vg_b_vol')
    E = f.read_variable('vg_e_vol')
    n = f.read_variable('proton/vg_rho')
    # Per components
    Bx = B[x1_global:x2_global,0,y1_global:y2_global,0]
    By = B[x1_global:x2_global,0,y1_global:y2_global,1]
    Bz = B[x1_global:x2_global,0,y1_global:y2_global,2]
    Ex = E[x1_global:x2_global,0,y1_global:y2_global,0]
    Ey = E[x1_global:x2_global,0,y1_global:y2_global,1]
    Ez = E[x1_global:x2_global,0,y1_global:y2_global,2]
    n  = n[cellids.argsort()].reshape([ysize_global,xsize_global])
    n  = n[y1_global:y2_global,x1_global:x2_global]

    # K_PERP
        # Slice to 1D. Will take values over a line in the middle of the box.
    Bx1D_kperp = BxR[int(BxR.shape[0]/2),:]
    By1D_kperp = ByR[int(ByR.shape[0]/2),:]
    Bz1D_kperp = BzR[int(BzR.shape[0]/2),:]
    Ex1D_kperp = ExR[int(ExR.shape[0]/2),:]
    Ey1D_kperp = EyR[int(EyR.shape[0]/2),:]
    Ez1D_kperp = EzR[int(EzR.shape[0]/2),:]
    n1D_kperp  = nR[int(nR.shape[0]/2),:]
        # Hanning windowing to force the edge of box to be periodic. Comment if boundary conditions are periodic.
    Bx1D_kperp *= np.hanning(len(Bx1D_kperp))
    By1D_kperp *= np.hanning(len(By1D_kperp))
    Bz1D_kperp *= np.hanning(len(Bz1D_kperp))
    Ex1D_kperp *= np.hanning(len(Ex1D_kperp))
    Ey1D_kperp *= np.hanning(len(Ey1D_kperp))
    Ez1D_kperp *= np.hanning(len(Ez1D_kperp))
    n1D_kperp  *= np.hanning(len(n1D_kperp))    
    # K_PARA
        # Slice to 1D. Will take values over a line in the middle of the box.
    Bx1D_kpara = BxR[:,int(BxR.shape[1]/2)]
    By1D_kpara = ByR[:,int(ByR.shape[1]/2)]
    Bz1D_kpara = BzR[:,int(BzR.shape[1]/2)]
    Ex1D_kpara = ExR[:,int(ExR.shape[1]/2)]
    Ey1D_kpara = EyR[:,int(EyR.shape[1]/2)]
    Ez1D_kpara = EzR[:,int(EzR.shape[1]/2)]
    n1D_kpara  = nR[:,int(nR.shape[1]/2)]
        # Hanning windowing to force the edge of box to be periodic. Comment if boundary conditions are periodic.
    Bx1D_kpara *= np.hanning(len(Bx1D_kpara)) 
    By1D_kpara *= np.hanning(len(By1D_kpara))
    Bz1D_kpara *= np.hanning(len(Bz1D_kpara))
    Ex1D_kpara *= np.hanning(len(Ex1D_kpara))
    Ey1D_kpara *= np.hanning(len(Ey1D_kpara))
    Ez1D_kpara *= np.hanning(len(Ez1D_kpara))
    n1D_kpara  *= np.hanning(len(n1D_kpara))

    out = [Bx1D_kperp,By1D_kperp,Bz1D_kperp,Ex1D_kperp,Ey1D_kperp,Ez1D_kperp,n1D_kperp,Bx1D_kpara,By1D_kpara,Bz1D_kpara,Ex1D_kpara,Ey1D_kpara,Ez1D_kpara,n1D_kpara]

    return(out)


def get_quantities():

    f = pt.vlsvfile.VlsvReader(filelist[int(len(filelist)/2)])
    cellids = f.read_variable('CellID')
    print('For '+filelist[int(len(filelist)/2)]+':')

    #Speed
    v     = f.read_variable("vg_v")
    v     = v[cellids.argsort()].reshape([ysize_global,xsize_global,3])
    v     = v[y1_global:y2_global,x1_global:x2_global,:]
    v     = [np.mean(v[:,:,0]),np.mean(v[:,:,1]),np.mean(v[:,:,2])]
    vnorm = np.linalg.norm(v)
    print('Vnorm = '+str(vnorm)+' m/s')

    #Bfield
    B     = f.read_fsgrid_variable('fg_b')
    B     = B[x1_global:x2_global,0,y1_global:y2_global,:]
    B     = [np.mean(B[:,:,0]),np.mean(B[:,:,1]),np.mean(B[:,:,2])]
    Bnorm = np.linalg.norm(B)
    print('Bnorm = '+str(Bnorm)+' T')

    if angle_global < 0.0:
        signe = - 1.0
    else:
        signe = 1.0

    vpara = signe * np.dot(B,v)/Bnorm
    print('Vpara = '+str(vpara)+' m/s')
    vperp = signe * np.linalg.norm(np.cross(B,np.cross(B,v))/Bnorm**2)
    print('Vperp = '+str(vperp)+' m/s')

    #Density
    rho = f.read_variable('proton/vg_rho')
    rho = rho[cellids.argsort()].reshape([ysize_global,xsize_global])
    rho = rho[y1_global:y2_global,x1_global:x2_global]
    rho = np.mean(rho)
    print('n = '+str(rho)+' m^-3')

    #Pressure
    P     = f.read_variable('proton/vg_ptensor_diagonal')
    P     = 1.0/3.0 * np.ma.sum(np.ma.masked_invalid(P),axis=-1)
    P     = P[cellids.argsort()].reshape([ysize_global,xsize_global])
    P     = P[y1_global:y2_global,x1_global:x2_global]
    Pnorm = np.mean(P)
    print('Pnorm = '+str(Pnorm)+' Pa')

    #Constants (SI units)
    q        = 1.60217e-19
    mu0      = 1.256637061e-6
    c        = 299792458
    mp       = 1.6726219e-27
    kB       = 1.38064852e-23
    me       = 9.1093826e-31
    epsilon0 = 8.85418781762e-12
    gamma    = 5.0/3.0
    d_i      = np.sqrt(mp*epsilon0*c*c/(rho*q*q))
    print('d_i = '+str(d_i)+' m')
    w_ci     = q*Bnorm/mp
    print('w_ci = '+str(w_ci)+' s^-1')
    vA       = Bnorm / np.sqrt(mu0 * rho * mp)
    print('vA = '+str(vA)+' m/s')
    vS       = np.sqrt(gamma * Pnorm / (rho * mp))
    print('vS = '+str(vS)+' m/s')

    return (vpara,vperp,d_i,w_ci,vA,vS,Bnorm,mu0,rho,mp)


def plt_2DFFT(filedir=None,
               start=None, stop=None,
               outputdir=None,
               boxre=[],boxm=[],
               run=None,
               Bmin=None, Bmax=None,
               Emin=None, Emax=None,
               nmin=None, nmax=None,
               dispersion=None,
               dr_ref=None,
               omegamax=None,
               numproc=40,
               paper=None):

    ''' Plots a 2D FFT with dispersion curves.

    :kword filedir:     Provide directory where files are located
    :kword start:       Bulk number for starting the data analysis
    :kword stop:        Bulk number for ending the data analysis
    :kword outputdir:   Path to directory where output files are created.
                        If directory does not exist, it will be created. If the string does not end in a
                        forward slash, the final parti will be used as a perfix for the files.
    :kword boxm:        Zoom box extents [x0,x1,y0,y1] in metres (default and truncate to: whole simulation box)
    :kword boxre:       Zoom box extents [x0,x1,y0,y1] in Earth radii (default and truncate to: whole simulation box)
    :kword run:         run identifier, used for constructing output filename
    :kword dr_ref:      Simulation resolution for PLOTTING in case different of the simulation of the run used for data (in meters)
    :kword Bmin,Bmax:   min and max values for colour scale and colour bar. If no values are given,
                        min and max values for whole plot are used.
    :kword Emin,Emax:   min and max values for colour scale and colour bar. If no values are given,
                        min and max values for whole plot are used.
    :kword nmin,nmax:   min and max values for colour scale and colour bar. If no values are given,
                        min and max values for whole plot are used.
    :kword dispersion:  List of strings for dispersion relations (Alfven, FMW, SMW, bulkV)
    :kword omegamax:    Max frequency for plotting.
    :kword numproc:     Number of processors for parallelisation (default: 40)
    :kword paper:       'E' or 'B'. Produces 2 side by side plot of B or E in para and perp direction ready for publication.


    # Example usage:
    plot_2DFFT(filedir=fileLocation,
               start=1600,stop=2400,
               run='BCQ',
               boxre=[0.0,3.0,15.0,18.0],
               Bmin=1e-7,Bmax=1e-4,
               Emin=1e-4,Emax=1e-1,
               nmin=1e-11,nmax=1e-7,
               dr_ref=300000,
               dispersion=['Alfven','FMW','SMW','bulkV']) 

    '''

    global filelist,start_global,boxcoords_global,x1_global,x2_global,y1_global,y2_global,xsize_global,ysize_global,angle_global,plane_global

    # Checks if outpudir exists, creates it otherwise
    if outputdir==None:
        outputdir=os.path.expandvars('$HOME/Plots/')
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    # Input file list
    if ((filedir!=None) and (start!=None) and (stop!=None)):
        filelist = []
        for step in range(start,stop+1):
            filename = filedir+'bulk.'+str(step).rjust(7,'0')+'.vlsv'
            filelist.append(filename)
    else:
        print("ERROR: needs a bulk file directory and start/stop steps")
        return


    RE = 6371e3 # m
    # Select box
    if len(boxm)==4:
        boxcoords=boxm
    elif len(boxre)==4:
        boxcoords=[i*RE for i in boxre]
    else:
        print("ERROR: No box was specified")
        return

    f = pt.vlsvfile.VlsvReader(filelist[int(len(filelist)/2)])
       
        # Get size of simulation (in cells)
    xsize = f.read_parameter('xcells_ini')
    xmax = f.read_parameter("xmax")
    xmin = f.read_parameter("xmin")
    dx   = (xmax-xmin)/xsize

    zsize = f.read_parameter('zcells_ini')
    if zsize == 1:
        plane = 'equatorial'
        ysize = f.read_parameter('ycells_ini')
        ymax = f.read_parameter("ymax")
        ymin = f.read_parameter("ymin")
        dy   = (ymax-ymin)/ysize
    else:
        plane = 'polar'
        ysize = zsize
        ymax = f.read_parameter("zmax")
        ymin = f.read_parameter("zmin")
        dy   = (ymax-ymin)/ysize
    print('-----------------------------')
    print('Run: '+run+', plane: '+plane+', box = ['+str(boxcoords[0]/RE)+','+str(boxcoords[1]/RE)+','+str(boxcoords[2]/RE)+','+str(boxcoords[3]/RE)+'] RE, start = '+str(start)+', stop = '+str(stop))
    print('-----------------------------')

        # Box coordinates # m
    x1 = int((boxcoords[0] - xmin) / dx)
    x2 = int((boxcoords[1] - xmin) / dx)
    y1 = int((boxcoords[2] - ymin) / dy)
    y2 = int((boxcoords[3] - ymin) / dy)

    start_global    = start
    boxcoods_global = boxcoords
    x1_global       = x1
    x2_global       = x2
    y1_global       = y1
    y2_global       = y2
    xsize_global    = xsize
    ysize_global    = ysize
    plane_global    = plane

    # Computes angle between magnetic field and z-axis (or y-axis for equatorial plane)
    # ASSUME THAT THE ANGLE DOESNT VARY TOO MUCH DURING THE ANALYSIS
    if __name__ == 'plot_2DFFT':
        pool         = Pool(numproc)
        print('Computing angle ...')
        return_angle = pool.map(get_angle, range(start,stop+1))
    else:
        print("Problem with angle parallelization")
        return

    angle_avg = np.mean(return_angle)
    angle_std = np.std(return_angle)
    print('-----------------------------')
    print('Angle = '+str(angle_avg)+' +/- '+str(angle_std)+' degrees.')
    print('-----------------------------')
    angle_global = angle_avg

    # Get 1D arrays of magnetic field, electric field and density
    if __name__ == 'plot_2DFFT':
        pool          = Pool(numproc)
        print('Getting fields ...')
        return_fields = pool.map(get_fields, range(start,stop+1))
    else:
        print('Problem with fields parallelisation')
        return

    return_fields = np.array(return_fields)

    B_kperp = np.zeros([return_fields.shape[0],3,return_fields.shape[2]])
    E_kperp = np.zeros([return_fields.shape[0],3,return_fields.shape[2]])
    n_kperp = np.zeros([return_fields.shape[0],return_fields.shape[2]])
    
    B_kpara = np.zeros([return_fields.shape[0],3,return_fields.shape[2]])
    E_kpara = np.zeros([return_fields.shape[0],3,return_fields.shape[2]])
    n_kpara = np.zeros([return_fields.shape[0],return_fields.shape[2]])
    
    # Creating fields
    for i in range(0,return_fields.shape[0]):
        B_kperp[i,0,:] = return_fields[i,0,:]
        B_kperp[i,1,:] = return_fields[i,1,:]        
        B_kperp[i,2,:] = return_fields[i,2,:] 
        E_kperp[i,0,:] = return_fields[i,3,:]
        E_kperp[i,1,:] = return_fields[i,4,:]
        E_kperp[i,2,:] = return_fields[i,5,:]
        n_kperp[i,:]   = return_fields[i,6,:]
 
        B_kpara[i,0,:] = return_fields[i,7,:]
        B_kpara[i,1,:] = return_fields[i,8,:]
        B_kpara[i,2,:] = return_fields[i,9,:]
        E_kpara[i,0,:] = return_fields[i,10,:]
        E_kpara[i,1,:] = return_fields[i,11,:]
        E_kpara[i,2,:] = return_fields[i,12,:]
        n_kpara[i,:]   = return_fields[i,13,:]

    print('-----------------------------')
    print('Arrays built')
    print('-----------------------------')
    
    # 2D FFT
    k_omega_Bx_kpara = np.fft.fftshift(np.fft.fft2(B_kpara[:,0].T)) #Shift to be centered around 0. Transposed to have complex part.
    k_omega_By_kpara = np.fft.fftshift(np.fft.fft2(B_kpara[:,1].T))
    k_omega_Bz_kpara = np.fft.fftshift(np.fft.fft2(B_kpara[:,2].T))
    k_omega_Ex_kpara = np.fft.fftshift(np.fft.fft2(E_kpara[:,0].T))
    k_omega_Ey_kpara = np.fft.fftshift(np.fft.fft2(E_kpara[:,1].T))
    k_omega_Ez_kpara = np.fft.fftshift(np.fft.fft2(E_kpara[:,2].T))
    k_omega_n_kpara  = np.fft.fftshift(np.fft.fft2(n_kpara.T))
    k_omega_Bx_kperp = np.fft.fftshift(np.fft.fft2(B_kperp[:,0].T))
    k_omega_By_kperp = np.fft.fftshift(np.fft.fft2(B_kperp[:,1].T))
    k_omega_Bz_kperp = np.fft.fftshift(np.fft.fft2(B_kperp[:,2].T))
    k_omega_Ex_kperp = np.fft.fftshift(np.fft.fft2(E_kperp[:,0].T))
    k_omega_Ey_kperp = np.fft.fftshift(np.fft.fft2(E_kperp[:,1].T))
    k_omega_Ez_kperp = np.fft.fftshift(np.fft.fft2(E_kperp[:,2].T))
    k_omega_n_kperp  = np.fft.fftshift(np.fft.fft2(n_kperp.T))


    # Compute quantities based on the middle bulk file of list
    vpara,vperp,d_i,w_ci,vA,vS,Bnorm,mu0,rho,mp = get_quantities()

    
    # Plot fourier space
    dt       = 0.5                         # dt of simulation (sec)
    kmin     = -np.pi / dx * d_i 
    kmax     = np.pi / dx * d_i

    omegamin = 0
    if omegamax == None:
       omegamax = np.pi / dt / w_ci
    
    if dr_ref != None:
        kmin_plot = - np.pi / dr_ref * d_i  
        kmax_plot = np.pi / dr_ref * d_i
    else:
        kmin_plot = kmin
        kmax_plot = kmax

    #Dispersion relations
    def CFL(k): #CFL condition
        return dx/0.044304 * (k / d_i)   #Constante may change with simulation
    
    def AWaves_para(k): #Alfven waves
        return abs(Bnorm / np.sqrt(mu0 * rho * mp) * (abs(k)/d_i) + vpara*(k/d_i)) # - flow speed
    
    def AWaves_perp(k): #Alfven waves
        return abs( Bnorm / np.sqrt(mu0 * rho * mp) * (abs(k)/d_i) - vperp*(k/d_i))
    
    def FMW_kpara(k): #Fast Magnetosonic Waves
        return abs(np.sqrt( 1.0/2.0 *( vS**2 + vA**2 + np.sqrt( (vS**2 + vA**2)**2 - 4* vA**2 * vS**2 ))) * (abs(k)/d_i) + vpara*(k/d_i))
    
    def SMW_kpara(k): #Slow Magnetosonic Waves
        return abs(np.sqrt( 1.0/2.0 *( vS**2 + vA**2 - np.sqrt( (vS**2 + vA**2)**2 - 4* vA**2 * vS**2 ))) * (abs(k)/d_i) + vpara*(k/d_i))
    
    def Bulk_v_para(k): #Bulk Velocity
        return  abs(vpara*(k/d_i))
    
    def FMW_kperp(k): #Fast Magnetosonic Waves
        return  abs(np.sqrt( 1.0/2.0 *( vS**2 + vA**2 + np.sqrt( (vS**2 + vA**2)**2))) * (abs(k)/d_i) - vperp*(k/d_i))
    
    def SMW_kperp(k): #Slow Magnetosonic Waves
        return  abs(np.sqrt( 1.0/2.0 *( vS**2 + vA**2 - np.sqrt( (vS**2 + vA**2)**2))) * (abs(k)/d_i) - vperp*(k/d_i))
    
    def Bulk_v_perp(k): #Bulk Velocity
        return  abs(vperp*(k/d_i))

    a           = np.linspace(kmin, kmax, 2000)
    cfl         = list(map(abs,map(CFL,a)))/w_ci
    awaves_para = list(map(AWaves_para,a))/w_ci
    awaves_perp = list(map(AWaves_perp,a))/w_ci
    fmw_kperp   = list(map(FMW_kperp,a))/w_ci
    smw_kperp   = list(map(SMW_kperp,a))/w_ci
    fmw_kpara   = list(map(FMW_kpara,a))/w_ci
    smw_kpara   = list(map(SMW_kpara,a))/w_ci
    bv_para     = list(map(Bulk_v_para,a))/w_ci
    bv_perp     = list(map(Bulk_v_perp,a))/w_ci
    
    print('-----------------------------')
    print('Plotting ...')
    print('-----------------------------')

    if paper == None:
        title_save = '2D_FFT'
        # K_PARA
            # Bx
        plt.figure(figsize=[4.0*3.5,3.15*3.5],dpi=300) #Set figure size
        ax    = plt.subplot()
       
        if Bmin != None:
           image = ax.imshow(abs(k_omega_Bx_kpara[int((k_omega_Bx_kpara.shape[0])/2):,:]),vmin=Bmin,vmax=Bmax, norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto')
        else:
           image = ax.imshow(abs(k_omega_Bx_kpara[int((k_omega_Bx_kpara.shape[0])/2):,:]), norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto')
    
        plt.locator_params(axis='x', nbins=10)      #Set number of ticks
        plt.setp(ax.get_xticklabels(),fontsize=30) #Set fontsize of tick labels
        plt.setp(ax.get_yticklabels(),fontsize=30)
        
        plt.xlim(kmin_plot, kmax_plot)
        plt.ylim(omegamin, omegamax)
        plt.xlabel("$k_\\parallel * d_i$",fontsize='35')
        plt.ylabel("$\\frac{\\omega}{\\Omega_c}$",fontsize='35')
        
        plt.plot(a, cfl, linewidth=3, color='black',label='CFL condition')
        
        for i in range(0,len(dispersion)):
            if dispersion[i] == 'Alfven':
                plt.plot(a,awaves_para,linewidth=3,color='blue',label='Alfven Waves')
            elif dispersion[i] == 'FMW':    
                plt.plot(a,fmw_kpara,linewidth=3,color='red',linestyle='--',label='Fast MW')
            elif dispersion[i] == 'SMW':    
                plt.plot(a,smw_kpara,linewidth=3,color='red',linestyle=':',label='Slow MW')
            elif dispersion[i] == 'bulkV':    
                plt.plot(a,bv_para,linewidth=3,color='black',linestyle='--',label='V')
        
        plt.legend(bbox_to_anchor=(-0.1,1.02,1.4,0.12),loc='upper right',fontsize='15',ncol=len(dispersion)+1,mode='expand') #Legend
    
        #Colorbar
        cb=plt.colorbar(image)
        plt.setp(plt.getp(cb.ax.axes, 'yticklabels'),fontsize=30) #Set fontisze of colorbar ticks
        
        variable   = 'B_x'
        cb.set_label('$'+variable+'^2 (T^2)$',fontsize=30)
        plt.title('$'+variable+'$',fontsize='40')
        plt.savefig(outputdir+run+'_'+variable+'_'+title_save+'_kpara.png')
        print('Figure saved at: '+outputdir+run+'_'+variable+'_'+title_save+'__kpara.png')
    
            # By
        plt.figure(figsize=[4.0*3.5,3.15*3.5],dpi=300) #Set figure size
        ax    = plt.subplot()
       
        if Bmin != None:
            image = ax.imshow(abs(k_omega_By_kpara[int((k_omega_By_kpara.shape[0])/2):,:]),vmin=Bmin,vmax=Bmax, norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto')
        else:
            image = ax.imshow(abs(k_omega_By_kpara[int((k_omega_By_kpara.shape[0])/2):,:]), norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto')
    
        plt.locator_params(axis='x', nbins=10)      #Set number of ticks
        plt.setp(ax.get_xticklabels(),fontsize=30) #Set fontsize of tick labels
        plt.setp(ax.get_yticklabels(),fontsize=30)
    
        plt.xlim(kmin_plot, kmax_plot)
        plt.ylim(omegamin, omegamax)
        plt.xlabel("$k_\\parallel * d_i$",fontsize='35')
        plt.ylabel("$\\frac{\\omega}{\\Omega_c}$",fontsize='35')
    
        plt.plot(a, cfl, linewidth=3, color='black',label='CFL condition')
    
        for i in range(0,len(dispersion)):
            if dispersion[i] == 'Alfven':
                plt.plot(a,awaves_para,linewidth=3,color='blue',label='Alfven Waves')
            elif dispersion[i] == 'FMW':
                plt.plot(a,fmw_kpara,linewidth=3,color='red',linestyle='--',label='Fast MW')
            elif dispersion[i] == 'SMW':
                plt.plot(a,smw_kpara,linewidth=3,color='red',linestyle=':',label='Slow MW')
            elif dispersion[i] == 'bulkV':
                plt.plot(a,bv_para,linewidth=3,color='black',linestyle='--',label='V')
    
        plt.legend(bbox_to_anchor=(-0.1,1.02,1.4,0.12),loc='upper right',fontsize='15',ncol=len(dispersion)+1,mode='expand') #Legend
    
        #Colorbar
        cb=plt.colorbar(image)
        plt.setp(plt.getp(cb.ax.axes, 'yticklabels'),fontsize=30) #Set fontisze of colorbar ticks
    
        variable   = 'B_y'
        cb.set_label('$'+variable+'^2 (T^2)$',fontsize=30)
        plt.title('$'+variable+'$',fontsize='40')
        plt.savefig(outputdir+run+'_'+variable+'_'+title_save+'_kpara.png')
        print ('Figure saved at: '+outputdir+run+'_'+variable+'_'+title_save+'_kpara.png')
    
        # Bz
        plt.figure(figsize=[4.0*3.5,3.15*3.5],dpi=300) #Set figure size
        ax    = plt.subplot()
    
        if Bmin != None:
            image = ax.imshow(abs(k_omega_Bz_kpara[int((k_omega_Bz_kpara.shape[0])/2):,:]),vmin=Bmin,vmax=Bmax, norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto')
        else:
            image = ax.imshow(abs(k_omega_Bz_kpara[int((k_omega_Bz_kpara.shape[0])/2):,:]), norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto')
    
        plt.locator_params(axis='x', nbins=10)      #Set number of ticks
        plt.setp(ax.get_xticklabels(),fontsize=30) #Set fontsize of tick labels
        plt.setp(ax.get_yticklabels(),fontsize=30)
    
        plt.xlim(kmin_plot, kmax_plot)
        plt.ylim(omegamin, omegamax)
        plt.xlabel("$k_\\parallel * d_i$",fontsize='35')
        plt.ylabel("$\\frac{\\omega}{\\Omega_c}$",fontsize='35')
    
        plt.plot(a, cfl, linewidth=3, color='black',label='CFL condition')
    
        for i in range(0,len(dispersion)):
            if dispersion[i] == 'Alfven':
                plt.plot(a,awaves_para,linewidth=3,color='blue',label='Alfven Waves')
            elif dispersion[i] == 'FMW':
                plt.plot(a,fmw_kpara,linewidth=3,color='red',linestyle='--',label='Fast MW')
            elif dispersion[i] == 'SMW':
                plt.plot(a,smw_kpara,linewidth=3,color='red',linestyle=':',label='Slow MW')
            elif dispersion[i] == 'bulkV':
                plt.plot(a,bv_para,linewidth=3,color='black',linestyle='--',label='V')
    
        plt.legend(bbox_to_anchor=(-0.1,1.02,1.4,0.12),loc='upper right',fontsize='15',ncol=len(dispersion)+1,mode='expand') #Legend
    
        #Colorbar
        cb=plt.colorbar(image)
        plt.setp(plt.getp(cb.ax.axes, 'yticklabels'),fontsize=30) #Set fontisze of colorbar ticks
    
        variable   = 'B_z'
        cb.set_label('$'+variable+'^2 (T^2)$',fontsize=30)
        plt.title('$'+variable+'$',fontsize='40')
        plt.savefig(outputdir+run+'_'+variable+'_'+title_save+'_kpara.png')
        print ('Figure saved at: '+outputdir+run+'_'+variable+'_'+title_save+'_kpara.png')
    
        # Ex
        plt.figure(figsize=[4.0*3.5,3.15*3.5],dpi=300) #Set figure size
        ax    = plt.subplot()
    
        if Emin != None:
            image = ax.imshow(abs(k_omega_Ex_kpara[int((k_omega_Ex_kpara.shape[0])/2):,:]),vmin=Emin,vmax=Emax, norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto')
        else:
            image = ax.imshow(abs(k_omega_Ex_kpara[int((k_omega_Ex_kpara.shape[0])/2):,:]), norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto')
    
        plt.locator_params(axis='x', nbins=10)      #Set number of ticks
        plt.setp(ax.get_xticklabels(),fontsize=30) #Set fontsize of tick labels
        plt.setp(ax.get_yticklabels(),fontsize=30)
    
        plt.xlim(kmin_plot, kmax_plot)
        plt.ylim(omegamin, omegamax)
        plt.xlabel("$k_\\parallel * d_i$",fontsize='35')
        plt.ylabel("$\\frac{\\omega}{\\Omega_c}$",fontsize='35')
    
        plt.plot(a, cfl, linewidth=3, color='black',label='CFL condition')
    
        for i in range(0,len(dispersion)):
            if dispersion[i] == 'Alfven':
                plt.plot(a,awaves_para,linewidth=3,color='blue',label='Alfven Waves')
            elif dispersion[i] == 'FMW':
                plt.plot(a,fmw_kpara,linewidth=3,color='red',linestyle='--',label='Fast MW')
            elif dispersion[i] == 'SMW':
                plt.plot(a,smw_kpara,linewidth=3,color='red',linestyle=':',label='Slow MW')
            elif dispersion[i] == 'bulkV':
                plt.plot(a,bv_para,linewidth=3,color='black',linestyle='--',label='V')
    
        plt.legend(bbox_to_anchor=(-0.1,1.02,1.4,0.12),loc='upper right',fontsize='15',ncol=len(dispersion)+1,mode='expand') #Legend
    
        #Colorbar
        cb=plt.colorbar(image)
        plt.setp(plt.getp(cb.ax.axes, 'yticklabels'),fontsize=30) #Set fontisze of colorbar ticks
    
        variable   = 'E_x'
        cb.set_label('$'+variable+'^2 (V/m)^2$',fontsize=30)
        plt.title('$'+variable+'$',fontsize='40')
        plt.savefig(outputdir+run+'_'+variable+'_'+title_save+'_kpara.png')
        print ('Figure saved at: '+outputdir+run+'_'+variable+'_'+title_save+'_kpara.png')
    
        # Ey
        plt.figure(figsize=[4.0*3.5,3.15*3.5],dpi=300) #Set figure size
        ax    = plt.subplot()
    
        if Emin != None:
            image = ax.imshow(abs(k_omega_Ey_kpara[int((k_omega_Ey_kpara.shape[0])/2):,:]),vmin=Emin,vmax=Emax, norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto')
        else:
            image = ax.imshow(abs(k_omega_Ey_kpara[int((k_omega_Ey_kpara.shape[0])/2):,:]), norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto')
    
        plt.locator_params(axis='x', nbins=10)      #Set number of ticks
        plt.setp(ax.get_xticklabels(),fontsize=30) #Set fontsize of tick labels
        plt.setp(ax.get_yticklabels(),fontsize=30)
    
        plt.xlim(kmin_plot, kmax_plot)
        plt.ylim(omegamin, omegamax)
        plt.xlabel("$k_\\parallel * d_i$",fontsize='35')
        plt.ylabel("$\\frac{\\omega}{\\Omega_c}$",fontsize='35')
    
        plt.plot(a, cfl, linewidth=3, color='black',label='CFL condition')
    
        for i in range(0,len(dispersion)):
            if dispersion[i] == 'Alfven':
                plt.plot(a,awaves_para,linewidth=3,color='blue',label='Alfven Waves')
            elif dispersion[i] == 'FMW':
                plt.plot(a,fmw_kpara,linewidth=3,color='red',linestyle='--',label='Fast MW')
            elif dispersion[i] == 'SMW':
                plt.plot(a,smw_kpara,linewidth=3,color='red',linestyle=':',label='Slow MW')
            elif dispersion[i] == 'bulkV':
                plt.plot(a,bv_para,linewidth=3,color='black',linestyle='--',label='V')
    
        plt.legend(bbox_to_anchor=(-0.1,1.02,1.4,0.12),loc='upper right',fontsize='15',ncol=len(dispersion)+1,mode='expand') #Legend
    
        #Colorbar
        cb=plt.colorbar(image)
        plt.setp(plt.getp(cb.ax.axes, 'yticklabels'),fontsize=30) #Set fontisze of colorbar ticks
    
        variable   = 'E_y'
        cb.set_label('$'+variable+'^2 (V/m)^2$',fontsize=30)
        plt.title('$'+variable+'$',fontsize='40')
        plt.savefig(outputdir+run+'_'+variable+'_'+title_save+'_kpara.png')
        print ('Figure saved at: '+outputdir+run+'_'+variable+'_'+title_save+'_kpara.png')
    
        # Ez
        plt.figure(figsize=[4.0*3.5,3.15*3.5],dpi=300) #Set figure size
        ax    = plt.subplot()
    
        if Emin != None:
            image = ax.imshow(abs(k_omega_Ez_kpara[int((k_omega_Ez_kpara.shape[0])/2):,:]),vmin=Emin,vmax=Emax, norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto')
        else:
            image = ax.imshow(abs(k_omega_Ez_kpara[int((k_omega_Ez_kpara.shape[0])/2):,:]), norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto')
    
        plt.locator_params(axis='x', nbins=10)      #Set number of ticks
        plt.setp(ax.get_xticklabels(),fontsize=30) #Set fontsize of tick labels
        plt.setp(ax.get_yticklabels(),fontsize=30)
    
        plt.xlim(kmin_plot, kmax_plot)
        plt.ylim(omegamin, omegamax)
        plt.xlabel("$k_\\parallel * d_i$",fontsize='35')
        plt.ylabel("$\\frac{\\omega}{\\Omega_c}$",fontsize='35')
    
        plt.plot(a, cfl, linewidth=3, color='black',label='CFL condition')
    
        for i in range(0,len(dispersion)):
            if dispersion[i] == 'Alfven':
                plt.plot(a,awaves_para,linewidth=3,color='blue',label='Alfven Waves')
            elif dispersion[i] == 'FMW':
                plt.plot(a,fmw_kpara,linewidth=3,color='red',linestyle='--',label='Fast MW')
            elif dispersion[i] == 'SMW':
                plt.plot(a,smw_kpara,linewidth=3,color='red',linestyle=':',label='Slow MW')
            elif dispersion[i] == 'bulkV':
                plt.plot(a,bv_para,linewidth=3,color='black',linestyle='--',label='V')
    
        plt.legend(bbox_to_anchor=(-0.1,1.02,1.4,0.12),loc='upper right',fontsize='15',ncol=len(dispersion)+1,mode='expand') #Legend
    
        #Colorbar
        cb=plt.colorbar(image)
        plt.setp(plt.getp(cb.ax.axes, 'yticklabels'),fontsize=30) #Set fontisze of colorbar ticks
    
        variable   = 'E_z'
        cb.set_label('$'+variable+'^2 (V/m)^2$',fontsize=30)
        plt.title('$'+variable+'$',fontsize='40')
        plt.savefig(outputdir+run+'_'+variable+'_'+title_save+'_kpara.png')
        print ('Figure saved at: '+outputdir+run+'_'+variable+'_'+title_save+'_kpara.png')
    
        # n
        plt.figure(figsize=[4.0*3.5,3.15*3.5],dpi=300) #Set figure size
        ax    = plt.subplot()
    
        if nmin != None:
            image = ax.imshow(abs(k_omega_n_kpara[int((k_omega_n_kpara.shape[0])/2):,:]),vmin=nmin,vmax=nmax, norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto')
        else:
            image = ax.imshow(abs(k_omega_n_kpara[int((k_omega_n_kpara.shape[0])/2):,:]), norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto')
    
        plt.locator_params(axis='x', nbins=10)      #Set number of ticks
        plt.setp(ax.get_xticklabels(),fontsize=30) #Set fontsize of tick labels
        plt.setp(ax.get_yticklabels(),fontsize=30)
    
        plt.xlim(kmin_plot, kmax_plot)
        plt.ylim(omegamin, omegamax)
        plt.xlabel("$k_\\parallel * d_i$",fontsize='35')
        plt.ylabel("$\\frac{\\omega}{\\Omega_c}$",fontsize='35')
    
        plt.plot(a, cfl, linewidth=3, color='black',label='CFL condition')
    
        for i in range(0,len(dispersion)):
            if dispersion[i] == 'Alfven':
                plt.plot(a,awaves_para,linewidth=3,color='blue',label='Alfven Waves')
            elif dispersion[i] == 'FMW':
                plt.plot(a,fmw_kpara,linewidth=3,color='red',linestyle='--',label='Fast MW')
            elif dispersion[i] == 'SMW':
                plt.plot(a,smw_kpara,linewidth=3,color='red',linestyle=':',label='Slow MW')
            elif dispersion[i] == 'bulkV':
                plt.plot(a,bv_para,linewidth=3,color='black',linestyle='--',label='V')
    
        plt.legend(bbox_to_anchor=(-0.1,1.02,1.4,0.12),loc='upper right',fontsize='15',ncol=len(dispersion)+1,mode='expand') #Legend
    
        #Colorbar
        cb=plt.colorbar(image)
        plt.setp(plt.getp(cb.ax.axes, 'yticklabels'),fontsize=30) #Set fontisze of colorbar ticks
    
        variable   = 'n'
        cb.set_label('$'+variable+'^2 (m^{-6})$',fontsize=30)
        plt.title('$'+variable+'$',fontsize='40')
        plt.savefig(outputdir+run+'_'+variable+'_'+title_save+'_kpara.png')
        print ('Figure saved at: '+outputdir+run+'_'+variable+'_'+title_save+'_kpara.png')
    
    # ----------------------------------------------------------------------------------------------------------------------------------------
    
        # K_PERP
            # Bx
        plt.figure(figsize=[4.0*3.5,3.15*3.5],dpi=300) #Set figure size
        ax    = plt.subplot()
       
        if Bmin != None:
           image = ax.imshow(abs(k_omega_Bx_kperp[int((k_omega_Bx_kperp.shape[0])/2):,:]),vmin=Bmin,vmax=Bmax, norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto')
        else:
           image = ax.imshow(abs(k_omega_Bx_kperp[int((k_omega_Bx_kperp.shape[0])/2):,:]), norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto')
    
        plt.locator_params(axis='x', nbins=10)      #Set number of ticks
        plt.setp(ax.get_xticklabels(),fontsize=30) #Set fontsize of tick labels
        plt.setp(ax.get_yticklabels(),fontsize=30)
        
        plt.xlim(kmin_plot, kmax_plot)
        plt.ylim(omegamin, omegamax)
        plt.xlabel("$k_\\parallel * d_i$",fontsize='35')
        plt.ylabel("$\\frac{\\omega}{\\Omega_c}$",fontsize='35')
        
        plt.plot(a, cfl, linewidth=3, color='black',label='CFL condition')
        
        for i in range(0,len(dispersion)):
            if dispersion[i] == 'Alfven':
                plt.plot(a,awaves_para,linewidth=3,color='blue',label='Alfven Waves')
            elif dispersion[i] == 'FMW':    
                plt.plot(a,fmw_kperp,linewidth=3,color='red',linestyle='--',label='Fast MW')
            elif dispersion[i] == 'SMW':    
                plt.plot(a,smw_kperp,linewidth=3,color='red',linestyle=':',label='Slow MW')
            elif dispersion[i] == 'bulkV':    
                plt.plot(a,bv_para,linewidth=3,color='black',linestyle='--',label='V')
        
        plt.legend(bbox_to_anchor=(-0.1,1.02,1.4,0.12),loc='upper right',fontsize='15',ncol=len(dispersion)+1,mode='expand') #Legend
    
        #Colorbar
        cb=plt.colorbar(image)
        plt.setp(plt.getp(cb.ax.axes, 'yticklabels'),fontsize=30) #Set fontisze of colorbar ticks
        
        variable   = 'B_x'
        cb.set_label('$'+variable+'^2 (T^2)$',fontsize=30)
        plt.title('$'+variable+'$',fontsize='40')
        plt.savefig(outputdir+run+'_'+variable+'_'+title_save+'_kperp.png')
        print ('Figure saved at: '+outputdir+run+'_'+variable+'_'+title_save+'_kperp.png')
    
            # By
        plt.figure(figsize=[4.0*3.5,3.15*3.5],dpi=300) #Set figure size
        ax    = plt.subplot()
       
        if Bmin != None:
            image = ax.imshow(abs(k_omega_By_kperp[int((k_omega_By_kperp.shape[0])/2):,:]),vmin=Bmin,vmax=Bmax, norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto')
        else:
            image = ax.imshow(abs(k_omega_By_kperp[int((k_omega_By_kperp.shape[0])/2):,:]), norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto')
    
        plt.locator_params(axis='x', nbins=10)      #Set number of ticks
        plt.setp(ax.get_xticklabels(),fontsize=30) #Set fontsize of tick labels
        plt.setp(ax.get_yticklabels(),fontsize=30)
    
        plt.xlim(kmin_plot, kmax_plot)
        plt.ylim(omegamin, omegamax)
        plt.xlabel("$k_\\parallel * d_i$",fontsize='35')
        plt.ylabel("$\\frac{\\omega}{\\Omega_c}$",fontsize='35')
    
        plt.plot(a, cfl, linewidth=3, color='black',label='CFL condition')
    
        for i in range(0,len(dispersion)):
            if dispersion[i] == 'Alfven':
                plt.plot(a,awaves_para,linewidth=3,color='blue',label='Alfven Waves')
            elif dispersion[i] == 'FMW':
                plt.plot(a,fmw_kperp,linewidth=3,color='red',linestyle='--',label='Fast MW')
            elif dispersion[i] == 'SMW':
                plt.plot(a,smw_kperp,linewidth=3,color='red',linestyle=':',label='Slow MW')
            elif dispersion[i] == 'bulkV':
                plt.plot(a,bv_para,linewidth=3,color='black',linestyle='--',label='V')
    
        plt.legend(bbox_to_anchor=(-0.1,1.02,1.4,0.12),loc='upper right',fontsize='15',ncol=len(dispersion)+1,mode='expand') #Legend
    
        #Colorbar
        cb=plt.colorbar(image)
        plt.setp(plt.getp(cb.ax.axes, 'yticklabels'),fontsize=30) #Set fontisze of colorbar ticks
    
        variable   = 'B_y'
        cb.set_label('$'+variable+'^2 (T^2)$',fontsize=30)
        plt.title('$'+variable+'$',fontsize='40')
        plt.savefig(outputdir+run+'_'+variable+'_'+title_save+'_kperp.png')
        print ('Figure saved at: '+outputdir+run+'_'+variable+'_'+title_save+'_kperp.png')
    
        # Bz
        plt.figure(figsize=[4.0*3.5,3.15*3.5],dpi=300) #Set figure size
        ax    = plt.subplot()
    
        if Bmin != None:
            image = ax.imshow(abs(k_omega_Bz_kperp[int((k_omega_Bz_kperp.shape[0])/2):,:]),vmin=Bmin,vmax=Bmax, norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto')
        else:
            image = ax.imshow(abs(k_omega_Bz_kperp[int((k_omega_Bz_kperp.shape[0])/2):,:]), norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto')
    
        plt.locator_params(axis='x', nbins=10)      #Set number of ticks
        plt.setp(ax.get_xticklabels(),fontsize=30) #Set fontsize of tick labels
        plt.setp(ax.get_yticklabels(),fontsize=30)
    
        plt.xlim(kmin_plot, kmax_plot)
        plt.ylim(omegamin, omegamax)
        plt.xlabel("$k_\\parallel * d_i$",fontsize='35')
        plt.ylabel("$\\frac{\\omega}{\\Omega_c}$",fontsize='35')
    
        plt.plot(a, cfl, linewidth=3, color='black',label='CFL condition')
    
        for i in range(0,len(dispersion)):
            if dispersion[i] == 'Alfven':
                plt.plot(a,awaves_para,linewidth=3,color='blue',label='Alfven Waves')
            elif dispersion[i] == 'FMW':
                plt.plot(a,fmw_kperp,linewidth=3,color='red',linestyle='--',label='Fast MW')
            elif dispersion[i] == 'SMW':
                plt.plot(a,smw_kperp,linewidth=3,color='red',linestyle=':',label='Slow MW')
            elif dispersion[i] == 'bulkV':
                plt.plot(a,bv_para,linewidth=3,color='black',linestyle='--',label='V')
    
        plt.legend(bbox_to_anchor=(-0.1,1.02,1.4,0.12),loc='upper right',fontsize='15',ncol=len(dispersion)+1,mode='expand') #Legend
    
        #Colorbar
        cb=plt.colorbar(image)
        plt.setp(plt.getp(cb.ax.axes, 'yticklabels'),fontsize=30) #Set fontisze of colorbar ticks
    
        variable   = 'B_z'
        cb.set_label('$'+variable+'^2 (T^2)$',fontsize=30)
        plt.title('$'+variable+'$',fontsize='40')
        plt.savefig(outputdir+run+'_'+variable+'_'+title_save+'_kperp.png')
        print ('Figure saved at: '+outputdir+run+'_'+variable+'_'+title_save+'_kperp.png')
    
        # Ex
        plt.figure(figsize=[4.0*3.5,3.15*3.5],dpi=300) #Set figure size
        ax    = plt.subplot()
    
        if Emin != None:
            image = ax.imshow(abs(k_omega_Ex_kperp[int((k_omega_Ex_kperp.shape[0])/2):,:]),vmin=Emin,vmax=Emax, norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto')
        else:
            image = ax.imshow(abs(k_omega_Ex_kperp[int((k_omega_Ex_kperp.shape[0])/2):,:]), norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto')
    
        plt.locator_params(axis='x', nbins=10)      #Set number of ticks
        plt.setp(ax.get_xticklabels(),fontsize=30) #Set fontsize of tick labels
        plt.setp(ax.get_yticklabels(),fontsize=30)
    
        plt.xlim(kmin_plot, kmax_plot)
        plt.ylim(omegamin, omegamax)
        plt.xlabel("$k_\\parallel * d_i$",fontsize='35')
        plt.ylabel("$\\frac{\\omega}{\\Omega_c}$",fontsize='35')
    
        plt.plot(a, cfl, linewidth=3, color='black',label='CFL condition')
    
        for i in range(0,len(dispersion)):
            if dispersion[i] == 'Alfven':
                plt.plot(a,awaves_para,linewidth=3,color='blue',label='Alfven Waves')
            elif dispersion[i] == 'FMW':
                plt.plot(a,fmw_kperp,linewidth=3,color='red',linestyle='--',label='Fast MW')
            elif dispersion[i] == 'SMW':
                plt.plot(a,smw_kperp,linewidth=3,color='red',linestyle=':',label='Slow MW')
            elif dispersion[i] == 'bulkV':
                plt.plot(a,bv_para,linewidth=3,color='black',linestyle='--',label='V')
    
        plt.legend(bbox_to_anchor=(-0.1,1.02,1.4,0.12),loc='upper right',fontsize='15',ncol=len(dispersion)+1,mode='expand') #Legend
    
        #Colorbar
        cb=plt.colorbar(image)
        plt.setp(plt.getp(cb.ax.axes, 'yticklabels'),fontsize=30) #Set fontisze of colorbar ticks
    
        variable   = 'E_x'
        cb.set_label('$'+variable+'^2 (V/m)^2$',fontsize=30)
        plt.title('$'+variable+'$',fontsize='40')
        plt.savefig(outputdir+run+'_'+variable+'_'+title_save+'_kperp.png')
        print ('Figure saved at: '+outputdir+run+'_'+variable+'_'+title_save+'_kperp.png')
    
        # Ey
        plt.figure(figsize=[4.0*3.5,3.15*3.5],dpi=300) #Set figure size
        ax    = plt.subplot()
    
        if Emin != None:
            image = ax.imshow(abs(k_omega_Ey_kperp[int((k_omega_Ey_kperp.shape[0])/2):,:]),vmin=Emin,vmax=Emax, norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto')
        else:
            image = ax.imshow(abs(k_omega_Ey_kperp[int((k_omega_Ey_kperp.shape[0])/2):,:]), norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto')
    
        plt.locator_params(axis='x', nbins=10)      #Set number of ticks
        plt.setp(ax.get_xticklabels(),fontsize=30) #Set fontsize of tick labels
        plt.setp(ax.get_yticklabels(),fontsize=30)
    
        plt.xlim(kmin_plot, kmax_plot)
        plt.ylim(omegamin, omegamax)
        plt.xlabel("$k_\\parallel * d_i$",fontsize='35')
        plt.ylabel("$\\frac{\\omega}{\\Omega_c}$",fontsize='35')
    
        plt.plot(a, cfl, linewidth=3, color='black',label='CFL condition')
    
        for i in range(0,len(dispersion)):
            if dispersion[i] == 'Alfven':
                plt.plot(a,awaves_para,linewidth=3,color='blue',label='Alfven Waves')
            elif dispersion[i] == 'FMW':
                plt.plot(a,fmw_kperp,linewidth=3,color='red',linestyle='--',label='Fast MW')
            elif dispersion[i] == 'SMW':
                plt.plot(a,smw_kperp,linewidth=3,color='red',linestyle=':',label='Slow MW')
            elif dispersion[i] == 'bulkV':
                plt.plot(a,bv_para,linewidth=3,color='black',linestyle='--',label='V')
    
        plt.legend(bbox_to_anchor=(-0.1,1.02,1.4,0.12),loc='upper right',fontsize='15',ncol=len(dispersion)+1,mode='expand') #Legend
    
        #Colorbar
        cb=plt.colorbar(image)
        plt.setp(plt.getp(cb.ax.axes, 'yticklabels'),fontsize=30) #Set fontisze of colorbar ticks
    
        variable   = 'E_y'
        cb.set_label('$'+variable+'^2 (V/m)^2$',fontsize=30)
        plt.title('$'+variable+'$',fontsize='40')
        plt.savefig(outputdir+run+'_'+variable+'_'+title_save+'_kperp.png')
        print ('Figure saved at: '+outputdir+run+'_'+variable+'_'+title_save+'_kperp.png')
    
        # Ez
        plt.figure(figsize=[4.0*3.5,3.15*3.5],dpi=300) #Set figure size
        ax    = plt.subplot()
    
        if Emin != None:
            image = ax.imshow(abs(k_omega_Ez_kperp[int((k_omega_Ez_kperp.shape[0])/2):,:]),vmin=Emin,vmax=Emax, norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto')
        else:
            image = ax.imshow(abs(k_omega_Ez_kperp[int((k_omega_Ez_kperp.shape[0])/2):,:]), norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto')
    
        plt.locator_params(axis='x', nbins=10)      #Set number of ticks
        plt.setp(ax.get_xticklabels(),fontsize=30) #Set fontsize of tick labels
        plt.setp(ax.get_yticklabels(),fontsize=30)
    
        plt.xlim(kmin_plot, kmax_plot)
        plt.ylim(omegamin, omegamax)
        plt.xlabel("$k_\\parallel * d_i$",fontsize='35')
        plt.ylabel("$\\frac{\\omega}{\\Omega_c}$",fontsize='35')
    
        plt.plot(a, cfl, linewidth=3, color='black',label='CFL condition')
    
        for i in range(0,len(dispersion)):
            if dispersion[i] == 'Alfven':
                plt.plot(a,awaves_para,linewidth=3,color='blue',label='Alfven Waves')
            elif dispersion[i] == 'FMW':
                plt.plot(a,fmw_kperp,linewidth=3,color='red',linestyle='--',label='Fast MW')
            elif dispersion[i] == 'SMW':
                plt.plot(a,smw_kperp,linewidth=3,color='red',linestyle=':',label='Slow MW')
            elif dispersion[i] == 'bulkV':
                plt.plot(a,bv_para,linewidth=3,color='black',linestyle='--',label='V')
    
        plt.legend(bbox_to_anchor=(-0.1,1.02,1.4,0.12),loc='upper right',fontsize='15',ncol=len(dispersion)+1,mode='expand') #Legend
    
        #Colorbar
        cb=plt.colorbar(image)
        plt.setp(plt.getp(cb.ax.axes, 'yticklabels'),fontsize=30) #Set fontisze of colorbar ticks
    
        variable   = 'E_z'
        cb.set_label('$'+variable+'^2 (V/m)^2$',fontsize=30)
        plt.title('$'+variable+'$',fontsize='40')
        plt.savefig(outputdir+run+'_'+variable+'_'+title_save+'_kperp.png')
        print ('Figure saved at: '+outputdir+run+'_'+variable+'_'+title_save+'_kperp.png')
    
        # n
        plt.figure(figsize=[4.0*3.5,3.15*3.5],dpi=300) #Set figure size
        ax    = plt.subplot()
    
        if nmin != None:
            image = ax.imshow(abs(k_omega_n_kperp[int((k_omega_n_kperp.shape[0])/2):,:]),vmin=nmin,vmax=nmax, norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto')
        else:
            image = ax.imshow(abs(k_omega_n_kperp[int((k_omega_n_kperp.shape[0])/2):,:]), norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto')
    
        plt.locator_params(axis='x', nbins=10)      #Set number of ticks
        plt.setp(ax.get_xticklabels(),fontsize=30) #Set fontsize of tick labels
        plt.setp(ax.get_yticklabels(),fontsize=30)
    
        plt.xlim(kmin_plot, kmax_plot)
        plt.ylim(omegamin, omegamax)
        plt.xlabel("$k_\\parallel * d_i$",fontsize='35')
        plt.ylabel("$\\frac{\\omega}{\\Omega_c}$",fontsize='35')
    
        plt.plot(a, cfl, linewidth=3, color='black',label='CFL condition')
    
        for i in range(0,len(dispersion)):
            if dispersion[i] == 'Alfven':
                plt.plot(a,awaves_para,linewidth=3,color='blue',label='Alfven Waves')
            elif dispersion[i] == 'FMW':
                plt.plot(a,fmw_kperp,linewidth=3,color='red',linestyle='--',label='Fast MW')
            elif dispersion[i] == 'SMW':
                plt.plot(a,smw_kperp,linewidth=3,color='red',linestyle=':',label='Slow MW')
            elif dispersion[i] == 'bulkV':
                plt.plot(a,bv_para,linewidth=3,color='black',linestyle='--',label='V')
    
        plt.legend(bbox_to_anchor=(-0.1,1.02,1.4,0.12),loc='upper right',fontsize='15',ncol=len(dispersion)+1,mode='expand') #Legend
    
        #Colorbar
        cb=plt.colorbar(image)
        plt.setp(plt.getp(cb.ax.axes, 'yticklabels'),fontsize=30) #Set fontisze of colorbar ticks
    
        variable   = 'n'
        cb.set_label('$'+variable+'^2 (m^{-6})$',fontsize=30)
        plt.title('$'+variable+'$',fontsize='40')
        plt.savefig(outputdir+run+'_'+variable+'_'+title_save+'_kperp.png')
        print ('Figure saved at: '+outputdir+run+'_'+variable+'_'+title_save+'_kperp.png')

    elif paper == 'E':
        fig, (ax_para, ax_perp) = plt.subplots(1,2,sharey=True)
        labelsize = 11
        fontsize  = 15
        
        cmap = 'jet'
        #K_PARA
        image = ax_para.imshow(abs(k_omega_Ey_kpara[int((k_omega_Ey_kpara.shape[0])/2):,:]),vmin=Emin,vmax=Emax, norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto', cmap = plt.get_cmap(cmap))
        
        ax_para.locator_params(axis='x', nbins=5)      #Set number of ticks
        ax_para.tick_params(axis='both',labelsize=labelsize)   #Set fontsize of tick labels
        ax_para.xaxis.set_major_locator(MultipleLocator(0.5))
        
        ax_para.set_xlim(kmin_plot, kmax_plot)
        ax_para.set_ylim(omegamin, omegamax)
        ax_para.set_xlabel("$k_\\parallel * d_i$",fontsize=fontsize)
        ax_para.set_ylabel("$\\frac{\\omega}{\\Omega_c}$",fontsize=fontsize)
        
        ax_para.plot(a, cfl, linewidth=2, color='black',label='CFL condition')
        ax_para.plot(a,awaves_para,linewidth=2,color='blue',label='Alfven Waves')
        ax_para.plot(a,fmw_kpara,linewidth=2,color='red',linestyle='--',label='Fast MW')
        ax_para.plot(a,smw_kpara,linewidth=2,color='red',linestyle=':',label='Slow MW')
        ax_para.plot(a,bv_para,linewidth=2,color='black',linestyle='--',label='V')
        
        ax_para.text(kmax_plot - 0.3,2.7,'(c)',fontsize=12)
        
        ax_para.legend(bbox_to_anchor = [-0.03, 1.04, 2.105, 0.05], loc='upper right',ncol=5, mode='expand', fontsize=8.5)
        
        #K_PERP
        image = ax_perp.imshow(abs(k_omega_Ey_kperp[int((k_omega_Ey_kperp.shape[0])/2):,:]),vmin=Emin,vmax=Emax, norm=colors.LogNorm(), extent=[kmin, kmax, omegamax, omegamin], aspect='auto',cmap = plt.get_cmap(cmap))
        
        ax_perp.locator_params(axis='x', nbins=5) #Set number of ticks
        ax_perp.tick_params(axis='x',labelsize=labelsize) #Set fontsize of tick labels
        ax_perp.xaxis.set_major_locator(MultipleLocator(0.5))
        
        ax_perp.set_xlim(kmin_plot, kmax_plot)
        ax_perp.set_ylim(omegamin, omegamax)
        ax_perp.set_xlabel("$k_\\bot * d_i$",fontsize=fontsize-1)
        
        ax_perp.plot(a, cfl, linewidth=2, color='black',label='CFL condition')
        ax_perp.plot(a,awaves_perp,linewidth=2,color='blue',label='Alfven Speed')
        ax_perp.plot(a,fmw_kperp,linewidth=2,color='red',linestyle='--',label='Fast MW')
        ax_perp.plot(a,smw_kperp,linewidth=2,color='red',linestyle=':',label='Slow MW')
        ax_perp.plot(a,bv_perp,linewidth=2,color='black',linestyle='--',label='V')
        
        ax_perp.text(kmax_plot - 0.3,2.7,'(d)',fontsize=12)
        
        fig.suptitle('$\\Delta r = '+str(int(dx/1000))+'$ km',fontsize=16)
        
        #plt.tight_layout(rect=[-0.01, 0.0, 1, 0.92])
        fig.subplots_adjust(wspace=0.05)

        #Colorbar
        fig.subplots_adjust(right=0.85)
        cbar_ax = fig.add_axes([0.855, 0.11, 0.03, 0.77])
        fig.colorbar(image, cax=cbar_ax)
        cbar_ax.set_ylabel('$E_y^2 (V/m)^2$',fontsize=fontsize-2)
        
        
        title = 'FFT_'+run
        #plt.savefig(path_fig+title+'.eps',dpi=400)
        plt.savefig(outputdir+title+'.pdf',dpi=400)
        print(outputdir+title+'.pdf')

    else:
        print('ERROR: kword paper specified with no valid argument')
        return

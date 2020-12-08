import pytools as pt
import numpy as np
from multiprocessing import Pool
import os,sys


numproc = 20 

run = sys.argv[1]
bulkStart  = int(sys.argv[2])
bulkEnd    = int(sys.argv[3])

path_bulk = '/wrk/users/dubart/300_test/'+run+'/bulk/'
path_save = '/wrk/users/dubart/analysis/subgrid/300_test/data/'

RE = 6371e3 # m

loc = ['mid','right','left','top_mid','bot_mid','top_r','bot_r','top_l','bot_l']

if run == '300_T2_B1' or run == '300_T3_B1' or run == '300_T4_B1' or run == '300_T5_B1':

    dt = 5.0 # How many seconds to calculate diffusion coeff. Min = 0.5s. 2.0 * Gyroperiod

elif run == '300_T2_B2' or run == '300_T3_B2' or run == '300_T4_B2' or run == '300_T5_B2':

    dt = 7.0

elif run == '300_T2_B3' or run == '300_T3_B3' or run == '300_T4_B3' or run == '300_T5_B3':

    dt = 9.0

nbins = 40

for i in range(0,len(loc)):

    CellID = int(np.load(path_save+'CellID_'+loc[i]+'.npy'))


    # Get f(mu,t)
    def get_fmu(step):
    
        pt.calculations.pitch_angles_box(vlsvReader=pt.vlsvfile.VlsvReader(path_bulk+'bulk.'+str(step).rjust(7,'0')+'.vlsv'),
                                             cellid = CellID, nbins = nbins,
                                             cosine = True,
                                             outputfile = path_save+run+'/PA/mu'+str(step).rjust(7,'0')+'_'+str(i))
    
        return
    
    pool = Pool(numproc)
    print('Getting fmu ...')
    pool.map(get_fmu,range(bulkStart,bulkEnd+1))
    pool.terminate()
    print('Saved "mu'+str(bulkStart).rjust(7,'0')+'" to "mu'+str(bulkEnd).rjust(7,'0')+'" at '+path_save)
    
    
    # Get first and second order derivatives in mu
    def get_dfdmu_dfdmumu(step):
    
        lines = []
    
        filePA = open(path_save+run+'/PA/mu'+str(step).rjust(7,'0')+'_'+str(i),'r')
    
        for line in filePA:
            lines.append(line)
    
        angles   = np.fromstring(lines[3],dtype=float,sep=' ')
        binvalue = np.fromstring(lines[5],dtype=float,sep=' ')
    
        filePA.close()
    
        dmu = abs(angles[0]-angles[1])
    
        dfdmu       = np.zeros(len(binvalue))
        dfdmu[1:-1] = (binvalue[2:] - binvalue[0:-2]) / (2*dmu) # center
        dfdmu[0]    = (binvalue[1]  - binvalue[0])    / dmu     # fw for left edge
        dfdmu[-1]   = (binvalue[-1] - binvalue[-2])   / dmu     # bw for right edge
    
        dfdmumu       = np.zeros(len(binvalue))
        dfdmumu[1:-1] = (binvalue[2:] - 2.0*binvalue[1:-1] + binvalue[0:-2]) / (dmu*dmu)
        dfdmumu[0]    = (binvalue[2]  - 2.0*binvalue[1]    + binvalue[0])    / (dmu*dmu)
        dfdmumu[-1]   = (binvalue[-1] - 2.0*binvalue[-2]   + binvalue[-3])   / (dmu*dmu)
    
        return [dfdmu,dfdmumu]
    
    pool = Pool(numproc)
    print('Getting dfdmu and dfdmumu ...')
    data = pool.map(get_dfdmu_dfdmumu,range(bulkStart,bulkEnd+1))
    pool.terminate()
    data = np.array(data)
    
    dfdmu   = data[:,0,:]
    dfdmumu = data[:,1,:]
    
    np.save(path_save+run+'/dfdmu_'+str(i)+'.npy',dfdmu)
    np.save(path_save+run+'/dfdmumu_'+str(i)+'.npy',dfdmumu)
    
    # Get time derivative
    def get_dfdt(step):
    
        if step == bulkStart:

            # t = 0 
            lines = []

            filePA_prev = open(path_save+run+'/PA/mu'+str(step).rjust(7,'0')+'_'+str(i),'r')

            for line in filePA_prev:
                lines.append(line)

            angles_prev   = np.fromstring(lines[3],dtype=float,sep=' ')
            binvalue_prev = np.fromstring(lines[5],dtype=float,sep=' ')

            filePA_prev.close()


            # from t = 1 to t = 2*dt

            binvalue_tmp = np.zeros(nbins)
            for j in range(step+1,step+int(dt*2.0)):

                lines = []
    
                filePA = open(path_save+run+'/PA/mu'+str(j).rjust(7,'0')+'_'+str(i),'r')
    
                for line in filePA:
                    lines.append(line)
    
                angles   = np.fromstring(lines[3],dtype=float,sep=' ')
                binvalue = np.fromstring(lines[5],dtype=float,sep=' ')
    
                filePA.close()

                binvalue_tmp = binvalue_tmp + binvalue

            binvalue_next = binvalue_tmp / int(dt*2.0)
 
            dfdt = (binvalue_next - binvalue_prev) / (dt/2.0)


        elif step > bulkStart and step < (bulkStart+int(dt*2.0)):

            binvalue_tmp = np.zeros(nbins)
            for j in range(bulkStart,step):

                lines = []

                filePA = open(path_save+run+'/PA/mu'+str(j).rjust(7,'0')+'_'+str(i),'r')

                for line in filePA:
                    lines.append(line)

                angles   = np.fromstring(lines[3],dtype=float,sep=' ')
                binvalue = np.fromstring(lines[5],dtype=float,sep=' ')

                filePA.close()

                binvalue_tmp = binvalue_tmp + binvalue

            binvalue_prev = binvalue_tmp / (step - bulkStart)

            binvalue_tmp = np.zeros(nbins)
            for j in range(step,(step+int(dt*2.0))):

                lines = []

                filePA = open(path_save+run+'/PA/mu'+str(j).rjust(7,'0')+'_'+str(i),'r')

                for line in filePA:
                    lines.append(line)

                angles   = np.fromstring(lines[3],dtype=float,sep=' ')
                binvalue = np.fromstring(lines[5],dtype=float,sep=' ')

                filePA.close()

                binvalue_tmp = binvalue_tmp + binvalue

            binvalue_next = binvalue_tmp / int(dt*2.0)

            dfdt = (binvalue_next - binvalue_prev) / ((dt/2.0) + (0.5 * (step - bulkStart) / 2.0))

        
        elif step > (bulkEnd-int(dt*2.0)) and step < bulkEnd:

            binvalue_tmp = np.zeros(nbins)
            for j in range((step-int(dt*2.0)),step):

                lines = []

                filePA = open(path_save+run+'/PA/mu'+str(j).rjust(7,'0')+'_'+str(i),'r')

                for line in filePA:
                    lines.append(line)

                angles   = np.fromstring(lines[3],dtype=float,sep=' ')
                binvalue = np.fromstring(lines[5],dtype=float,sep=' ')

                filePA.close()

                binvalue_tmp = binvalue_tmp + binvalue

            binvalue_prev = binvalue_tmp / int(dt*2.0)

            binvalue_tmp = np.zeros(nbins)
            for j in range(step,bulkEnd):

                lines = []

                filePA = open(path_save+run+'/PA/mu'+str(j).rjust(7,'0')+'_'+str(i),'r')

                for line in filePA:
                    lines.append(line)

                angles   = np.fromstring(lines[3],dtype=float,sep=' ')
                binvalue = np.fromstring(lines[5],dtype=float,sep=' ')

                filePA.close()

                binvalue_tmp = binvalue_tmp + binvalue

            binvalue_next = binvalue_tmp / (bulkEnd - step)

            dfdt = (binvalue_next - binvalue_prev) / ((dt/2.0) + (0.5 * (bulkEnd - step) / 2.0))


        else:

            binvalue_tmp = np.zeros(nbins)
            for j in range((step-int(dt*2.0)),step):

                lines = []

                filePA = open(path_save+run+'/PA/mu'+str(j).rjust(7,'0')+'_'+str(i),'r')

                for line in filePA:
                    lines.append(line)

                angles   = np.fromstring(lines[3],dtype=float,sep=' ')
                binvalue = np.fromstring(lines[5],dtype=float,sep=' ')

                filePA.close()

                binvalue_tmp = binvalue_tmp + binvalue

            binvalue_prev = binvalue_tmp / int(dt*2.0)

            binvalue_tmp = np.zeros(nbins)
            for j in range(step,(step+int(dt*2.0))):

                lines = []

                filePA = open(path_save+run+'/PA/mu'+str(j).rjust(7,'0')+'_'+str(i),'r')

                for line in filePA:
                    lines.append(line)

                angles   = np.fromstring(lines[3],dtype=float,sep=' ')
                binvalue = np.fromstring(lines[5],dtype=float,sep=' ')

                filePA.close()

                binvalue_tmp = binvalue_tmp + binvalue

            binvalue_next = binvalue_tmp / int(dt*2.0)

            dfdt = (binvalue_next - binvalue_prev) / dt

    
        return [dfdt]
    
    pool = Pool(numproc)
    print('Getting dfdt ...')
    data = pool.map(get_dfdt,range(bulkStart,bulkEnd+1))
    pool.terminate()
    data = np.array(data)
    
    dfdt = data[:,0,:]
    
    np.save(path_save+run+'/dfdt_'+str(i)+'.npy',dfdt)
    
    # Get dmu
    lines = []
    
    filePA = open(path_save+run+'/PA/mu'+str(bulkStart).rjust(7,'0')+'_'+str(i),'r')
    
    for line in filePA:
        lines.append(line)
    
    angles   = np.fromstring(lines[3],dtype=float,sep=' ')
    
    filePA.close()
    
    dmu = abs(angles[0]-angles[1])
    
    # Get diffusion coefficient
    def get_Dmumu(step):
    
        Mat = np.zeros([dfdmumu.shape[1],dfdmumu.shape[1]])
    
        for j in range(0,dfdmumu.shape[1]):
            for k in range(0,dfdmumu.shape[1]):
    
                if k == j:
                    Mat[j,k] = dfdmumu[step,j] # Diagonal
                elif k == j+1:
                    Mat[j,k] = dfdmu[step,j]/(2*dmu) # upper
                elif k == j-1:
                    Mat[j,k] = -dfdmu[step,j]/(2*dmu) # lower
    
    
    #    Mat[0,-1] = -dfdmu[step,0]/(2*dmu) # added for periodic boundary conditions
    #    Mat[-1,0] = dfdmu[step,-1]/(2*dmu)
    
        Dmumu = np.linalg.solve(Mat,dfdt[step,:])
    
        return [Dmumu]
    
      
    pool = Pool(numproc)
    print('Getting Dmumu ...')
    data = pool.map(get_Dmumu,range(0,dfdt.shape[0]))
    pool.terminate()
    
    data = np.array(data)
    
    Dmumu = data[:,0,:]
    
    np.save(path_save+run+'/Dmumu_'+str(i)+'.npy',Dmumu)
    print('Saved '+path_save+run+'/Dmumu_'+str(i)+'.npy')

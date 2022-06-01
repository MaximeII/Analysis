import numpy as np
import pytools as pt
import matplotlib
import matplotlib.pyplot as plt
import sys
import scipy
from scipy.sparse.linalg import bicg

run = sys.argv[1]

if run == '300':
    path_bulk = '/wrk/users/dubart/dmumu/'+run+'/bulk/'
else:
    path_bulk = '/wrk/users/dubart/diff_test/dmumu/diff/'+run+'/CFL1.0/bulk/'
    #path_bulk = '/wrk/users/dubart/diff_test/dmumu/diff/'+run+'/opti_v1/bulk/'

#path_bulk = '/wrk/users/dubart/diff_test/'+run+'/mu_bulk/'
path_save = '/wrk/users/dubart/analysis/subgrid/300_test/Green_diff/data/mu/'+run+'/'

bulkStart = int(sys.argv[2])
bulkEnd   = int(sys.argv[3])
timetot   = bulkEnd - bulkStart +1

mode = sys.argv[4]

tsize = 4 #timetot

nbins = 30 # number of bins in mu

fmu_array = np.empty([nbins,tsize],dtype=np.float128)

for step in range(0,tsize):

    file_mu = open(path_save + 'mu_'+run+'_'+str(step*100).rjust(7,'0')+'.txt','r')
    Lines   = file_mu.readlines()
    file_mu.close()

    fmu = np.fromstring(Lines[5],dtype='float',sep=' ')

    fmu_array[:,step] = fmu

    if step == bulkStart:
        mu  = np.fromstring(Lines[3],dtype='float',sep=' ')
        mu_mid = (mu[:-1] + mu[1:])/2

dfdt    = np.empty([nbins,tsize],dtype=np.float128) # d
dfdmu   = np.empty([nbins,tsize],dtype=np.float128)
dfdmumu = np.empty([nbins,tsize],dtype=np.float128)

dmu = abs(mu_mid[1]- mu_mid[0]) 
dt  = 10

# dfdt
dfdt[:,0]    = (fmu_array[:,1] - fmu_array[:,0]) / dt         # Forward
dfdt[:,1:-1] = (fmu_array[:,2:] - fmu_array[:,:-2]) / (2*dt)  # Central
dfdt[:,-1]   = (fmu_array[:,-1] - fmu_array[:,-2]) / dt       # Backward

# dfdmu
dfdmu[0,:]  = (fmu_array[1,:] - fmu_array[0,:]) / dmu         
dfdmu[1:-1,:] = (fmu_array[2:,:] - fmu_array[:-2,:]) / (2*dmu)  
dfdmu[-1,:] = (fmu_array[-1,:] - fmu_array[-2,:]) / dmu

# dfdmumu
dfdmumu[0,:]    = (fmu_array[2,:] - 2*fmu_array[1,:] + fmu_array[0,:]) / dmu**2
dfdmumu[1:-1,:] = (fmu_array[2:,:] - 2*fmu_array[1:-1,:] + fmu_array[:-2,:]) / dmu**2
dfdmumu[-1,:]   = (fmu_array[-1,:] - 2*fmu_array[-2,:] + fmu_array[-3,:]) / dmu**2

Dmumu      = np.empty([nbins,tsize],dtype=np.float128)


indx_minus = np.where(mu <= -0.2)[0][-1]
indx_plus = np.where(mu >= 0.2)[0][0]

for step in range(0,tsize):

    try:
        if mode == 'half':

            # Matrix = np.diag(dfdmumu[:,step],0) + np.diag(dfdmu[:-1,step]/(2*dmu),1) + np.diag(-dfdmu[1:,step]/(2*dmu),-1)
            # 
            # InvM = scipy.linalg.inv(Matrix)

            # Dmumu[:,step] = np.matmul(InvM,dfdt[:,step])

            Matrix = np.diag(dfdmumu[:indx_minus+1,step],0) + np.diag(dfdmu[:indx_minus,step]/(2*dmu),1) + np.diag(-dfdmu[1:indx_minus+1,step]/(2*dmu),-1)

            #InvM = scipy.linalg.inv(Matrix)

            #Dmumu[:indx_minus+1,step] = np.matmul(InvM,dfdt[:indx_minus+1,step])
 
            Dmumu[:indx_minus+1,step] = scipy.linalg.solve(Matrix,dfdt[:indx_minus+1,step])

            Matrix = np.diag(dfdmumu[indx_plus:,step],0) + np.diag(dfdmu[indx_plus:-1,step]/(2*dmu),1) + np.diag(-dfdmu[indx_plus+1:,step]/(2*dmu),-1)

            #InvM = scipy.linalg.inv(Matrix)

            #Dmumu[indx_plus:,step] = np.matmul(InvM,dfdt[indx_plus:,step])

            Dmumu[indx_plus:,step] = scipy.linalg.solve(Matrix,dfdt[indx_plus:,step])

        elif mode == 'full':

            Matrix = np.diag(dfdmumu[:,step],0) + np.diag(dfdmu[:-1,step]/(2*dmu),1) + np.diag(-dfdmu[1:,step]/(2*dmu),-1)
            

            #InvM = scipy.linalg.inv(Matrix)

            #Dmumu[:,step] = np.matmul(InvM,dfdt[:,step])
            Dmumu[:,step] = scipy.linalg.solve(Matrix,dfdt[:,step])

            print(Dmumu[:,step])

    except Exception as e:
       print('i = '+str(step)+': '+str(e))

np.save(path_save+'Dmumu_'+run+'.npy',Dmumu)
print(path_save+'Dmumu_'+run+'.npy')

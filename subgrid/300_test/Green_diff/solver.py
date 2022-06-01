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

bulks = [bulkStart,bulkEnd]

nbins = 30 # number of bins in mu

fmu_array = np.empty([nbins,2],dtype=np.float128)

for step in range(0,len(bulks)):

    file_mu = open(path_save + 'mu_'+run+'_'+str(bulks[step]).rjust(7,'0')+'.txt','r')
    Lines   = file_mu.readlines()
    file_mu.close()

    fmu = np.fromstring(Lines[5],dtype='float',sep=' ')

    fmu_array[:,step] = fmu

    mu  = np.fromstring(Lines[3],dtype='float',sep=' ')
    mu_mid = (mu[:-1] + mu[1:])/2

fmu_avg = np.mean(fmu_array,axis=1)

dfdmu   = np.empty(nbins,dtype=np.float128)
dfdmumu = np.empty(nbins,dtype=np.float128)

dmu = abs(mu_mid[1]- mu_mid[0]) 
dt  = (timetot)*0.1

# dfdt
dfdt = (fmu_array[:,1] - fmu_array[:,0]) / dt

# dfdmu
dfdmu[0]    = (fmu_avg[1]  - fmu_avg[0])   / dmu         
dfdmu[1:-1] = (fmu_avg[2:] - fmu_avg[:-2]) / (2*dmu)  
dfdmu[-1]   = (fmu_avg[-1] - fmu_avg[-2])  / dmu

# dfdmumu
dfdmumu[0]    = (fmu_avg[2]  - 2*fmu_avg[1]    + fmu_avg[0])   / dmu**2
dfdmumu[1:-1] = (fmu_avg[2:] - 2*fmu_avg[1:-1] + fmu_avg[:-2]) / dmu**2
dfdmumu[-1]   = (fmu_avg[-1] - 2*fmu_avg[-2]   + fmu_avg[-3])  / dmu**2

try:
    Matrix = np.diag(dfdmumu,0) + np.diag(dfdmu[:-1]/(2*dmu),1) + np.diag(-dfdmu[1:]/(2*dmu),-1)
        
    Dmumu = scipy.linalg.solve(Matrix,dfdt)


except Exception as e:
   print('i = '+str(step)+': '+str(e))

np.save(path_save+'Dmumu_matrix_'+run+'_'+str(bulkStart)+'-'+str(bulkEnd)+'.npy',Dmumu)
print(path_save+'Dmumu_matrix_'+run+'_'+str(bulkStart)+'-'+str(bulkEnd)+'.npy')

dfdt_integration = np.empty(nbins,dtype=np.float128)

for i in range(0,len(dfdt)):
    dfdt_integration[i] = scipy.integrate.simps(dfdt[:i+1],x=mu_mid[:i+1],dx=dmu)

Dmumu = dfdt_integration/dfdmu

np.save(path_save+'Dmumu_integration_'+run+'_'+str(bulkStart)+'-'+str(bulkEnd)+'.npy',Dmumu)
print(path_save+'Dmumu_integration_'+run+'_'+str(bulkStart)+'-'+str(bulkEnd)+'.npy')


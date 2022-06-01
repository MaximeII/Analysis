import numpy as np
import pytools as pt
import matplotlib
import matplotlib.pyplot as plt
import sys

run = sys.argv[1]

if run == '300':
    path_bulk = '/wrk/users/dubart/dmumu/'+run+'/bulk/'
else:
    path_bulk = '/wrk/users/dubart/diff_test/dmumu/diff/'+run+'/CFL1.0/bulk/'
    #path_bulk = '/wrk/users/dubart/diff_test/dmumu/diff/'+run+'/opti_v1/bulk/'

path_save = '/wrk/users/dubart/analysis/subgrid/300_test/Green_diff/data/mu/'+run+'/'

bulkStart = int(sys.argv[2])
bulkEnd   = int(sys.argv[3])
timetot   = bulkEnd - bulkStart +1

nbins = 30 # number of bins in mu

#------------------------------------------------------------------
## Tri Diagonal Matrix Algorithm(a.k.a Thomas algorithm) solver
def TDMAsolver(a, b, c, d):
    '''
    TDMA solver, a b c d can be NumPy array type or Python list type.
    refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    and to http://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)
    '''
    nf = len(d) # number of equations
    ac, bc, cc, dc = map(np.array, (a, b, c, d)) # copy arrays
    for it in range(1, nf):
        mc = ac[it-1]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1] 
        dc[it] = dc[it] - mc*dc[it-1]
        	    
    xc = bc
    xc[-1] = dc[-1]/bc[-1]

    for il in range(nf-2, -1, -1):
        xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]

    return xc
#-----------------------------------------------------------------

fmu_array = np.empty([nbins,timetot])

for step in range(bulkStart,bulkEnd+1):

    file_mu = open(path_save + 'mu_'+run+'_'+str(step).rjust(7,'0')+'.txt','r')
    Lines   = file_mu.readlines()
    file_mu.close()

    fmu = np.fromstring(Lines[5],dtype='float',sep=' ')

    fmu_array[:,step] = fmu

    if step == bulkStart:
        mu  = np.fromstring(Lines[3],dtype='float',sep=' ')
        mu_mid = (mu[:-1] + mu[1:])/2  


dfdt    = np.empty([nbins,timetot]) # d
dfdmu   = np.empty([nbins,timetot])
dfdmumu = np.empty([nbins,timetot])

dmu = abs(mu_mid[1]- mu_mid[0]) 
dt  = 0.1

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

Dmumu      = np.empty([nbins,timetot])

for step in range(bulkStart,bulkEnd+1):

    Dmumu[:,step] = TDMAsolver(-dfdmu[1:,step]/(2*dmu),dfdmumu[:,step],dfdmu[:-1,step]/(2*dmu),dfdt[:,step])

np.save(path_save+'Dmumu_'+run+'.npy',Dmumu)
print(path_save+'Dmumu_'+run+'.npy')

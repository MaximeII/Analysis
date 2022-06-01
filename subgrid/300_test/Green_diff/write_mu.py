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

nbins = 30 # number of bins in mu

for step in range(bulkStart,bulkEnd+1):

    file_mu = open(path_save + 'mu_'+run+'_'+str(step).rjust(7,'0')+'.txt','r')
    Lines   = file_mu.readlines()
    file_mu.close()

    fmu = np.fromstring(Lines[5],dtype='float',sep=' ')

    mu  = np.fromstring(Lines[3],dtype='float',sep=' ')
    mu_mid = (mu[:-1] + mu[1:])/2

    data = np.column_stack([mu_mid,fmu])
    np.savetxt(path_save+'mu_column_'+str(step).rjust(7,'0')+'.txt',data)

    print(path_save+'mu_column_'+str(step).rjust(7,'0')+'.txt')


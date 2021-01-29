import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt
from multiprocessing import Pool

numproc = 20 

run = sys.argv[1]
bulkStart  = int(sys.argv[2])
bulkEnd    = int(sys.argv[3])

path_bulk = '/wrk/users/dubart/300_test/'+run+'/bulk/'
path_save = '/wrk/users/dubart/analysis/subgrid/300_test/data/'+run+'/'
path_fig  = '/wrk/users/dubart/analysis/subgrid/300_test/fig/'+run+'/PA_DC/'

loc = ['mid','right','left','top_mid','bot_mid','top_r','bot_r','top_l','bot_l']

RE = 6371e3 # m

# Initialise arrays
lines = []

filePA = open(path_save+'PA/mu'+str(bulkStart).rjust(7,'0')+'_0','r')

for line in filePA:
    lines.append(line)

angles   = np.fromstring(lines[3],dtype=float,sep=' ')
binvalue = np.fromstring(lines[5],dtype=float,sep=' ')

filePA.close()

angles_center = (angles[0:-1]+angles[1:])/2
length_bin    = len(binvalue)

def variance(step):

    binvalue_tmp = np.zeros(length_bin)
    
    for j in range(0,len(loc)):
    
        # f(mu,t)
        lines = []
    
        filePA = open(path_save+'PA/mu'+str(step).rjust(7,'0')+'_'+str(j),'r')
    
        for line in filePA:
            lines.append(line)
    
        binvalue = np.fromstring(lines[5],dtype=float,sep=' ')
    
        filePA.close()
    
        binvalue_tmp = binvalue_tmp + binvalue
    
    
    # Avg loc
    binvalue_avg = binvalue_tmp / len(loc)
    
    weighted_avg = np.average(angles_center,weights=binvalue_avg)
    bin_std      = np.sqrt(np.average((angles_center-weighted_avg)**2,weights=binvalue_avg))
    
    return bin_std

pool = Pool(numproc)
data = pool.map(variance,range(bulkStart,bulkEnd+1))
pool.terminate()
data = np.array(data)
np.save(path_save+'std_'+run+'.npy',data)

print(path_save+'std_'+run+'.npy')

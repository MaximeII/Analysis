import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt

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

diff_binvalue_tmp = np.zeros(length_bin)

for i in range(bulkStart,bulkEnd+1):

       binvalue_tmp_present = np.zeros(length_bin)
       binvalue_tmp_next    = np.zeros(length_bin)
       
       for j in range(0,len(loc)):

           # f(mu,t)
           lines = []

           filePA = open(path_save+'PA/mu'+str(i).rjust(7,'0')+'_'+str(j),'r')

           for line in filePA:
               lines.append(line)

           binvalue = np.fromstring(lines[5],dtype=float,sep=' ')

           filePA.close()

           binvalue_tmp_present = binvalue_tmp_present + binvalue
 
 
           # f(mu,t+1)
           lines = []

           filePA = open(path_save+'PA/mu'+str(i+1).rjust(7,'0')+'_'+str(j),'r')

           for line in filePA:
               lines.append(line)

           binvalue = np.fromstring(lines[5],dtype=float,sep=' ')

           filePA.close()

           binvalue_tmp_next = binvalue_tmp_next + binvalue

       # Avg loc
       binvalue_avg_present = binvalue_tmp_present / len(loc)
       binvalue_avg_next    = binvalue_tmp_next / len(loc)

       diff_binvalue_tmp = diff_binvalue_tmp + ((binvalue_avg_next - binvalue_avg_present) / ((binvalue_avg_present + binvalue_avg_next)/2.0))

diff_binvalue_avg = diff_binvalue_tmp / (bulkEnd - bulkStart + 1) 

np.save(path_save+'Diff_fmu_'+run+'.npy',diff_binvalue_avg)
print(path_save+'Diff_fmu_'+run+'.npy')

      
    

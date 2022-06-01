import pytools as pt
import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt

run = sys.argv[1]

if run == '300':
    path_bulk = '/wrk-vakka/users/dubart/diff_test/dmumu/'+run+'/bulk/'
elif run == 'Markus':
    path_bulk = '/wrk-vakka/users/markusb/weirddiffusion/Losscone_may31/'
elif run == 'Markus_artificial':
    path_bulk = '/wrk-vakka/users/markusb/weirddiffusion/Losscone_may31_artificial/'
else:
    path_bulk = '/wrk/users/dubart/diff_test/dmumu/diff/'+run+'/CFL1.0/bulk/'

#path_bulk = '/wrk/users/dubart/diff_test/dmumu/diff/'+run+'/opti_v1/bulk/'
#path_bulk = '/wrk-vakka/users/dubart/diff_test/'+run+'/mu_bulk/'
#path_bulk = '/wrk-vakka/users/dubart/diff_test/dmumu/'+run+'/bulk/'
path_save = '/wrk/users/dubart/analysis/subgrid/300_test/Green_diff/data/mu/'+run+'/'
path_fig  = '/wrk/users/dubart/analysis/subgrid/300_test/Green_diff/fig/mu/'+run+'/'


if len(sys.argv) == 4: # Starting and end frames given
    timetot = range(int(sys.argv[2]), int(sys.argv[3]),1)
else: # Only starting frame given, generate one frame
    timetot = range(int(sys.argv[2]), int(sys.argv[2])+1, 1)
for step in timetot:

    if step == 0:
        bulkname = "initial-grid.0000000.vlsv"
    else:
        # Source data file
        bulkname = "bulk."+str(step).rjust(7,'0')+".vlsv"
    print(bulkname)

    #CellID = int(np.load(path_save+'CellID_'+run+'.npy'))
    CellID = 1

    if run == 'Markus' or run == 'Markus_artificial':

       pt.calculations.pitch_angles(vlsvReader  = None,
                                    filename    = path_bulk + bulkname,
                                    cellid      = CellID,
                                    nbins       = 30,
                                    cosine      = True,
                                    plasmaframe = False,
                                    vcut        = 100e3,
                                    pop         = 'protonmono',
                                    outputfile  = path_save + 'mu_'+run+'_'+str(step).rjust(7,'0')+'.txt')

    else:

       pt.calculations.pitch_angles(vlsvReader  = None,
                                    filename    = path_bulk + bulkname,
                                    cellid      = CellID,
                                    nbins       = 30,
                                    cosine      = True,
                                    plasmaframe = True,
                                    vcut        = 100e3,
                                    #pop         = 'protonmono',
                                    outputfile  = path_save + 'mu_'+run+'_'+str(step).rjust(7,'0')+'.txt')
    print(path_save + 'mu_'+run+'_'+str(step).rjust(7,'0')+'.txt')

       

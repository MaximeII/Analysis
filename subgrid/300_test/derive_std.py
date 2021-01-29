import pytools as pt
import numpy as np
import sys

runs = ['300_T2_B1','300_T3_B1','300_T4_B1','300_T5_B1','300_T2_B2','300_T3_B2','300_T4_B2','300_T5_B2','300_T2_B3','300_T3_B3','300_T4_B3','300_T5_B3']

for i in range(0,len(runs)):

    path_save = '/wrk/users/dubart/analysis/subgrid/300_test/data/'+runs[i]+'/'

    std = np.load(path_save+'std_'+runs[i]+'.npy')

    dt = 0.5
    time = np.linspace(0.5,50,len(std))

    if runs[i][5] == '3':
        tmin = int(np.where(time == 20.0)[0])
        tmax = int(np.where(time == 30.0)[0])

    elif runs[i][5] == '4':
        tmin = int(np.where(time == 12.0)[0])
        tmax = int(np.where(time == 22.0)[0])

    elif runs[i][5] == '5':
        tmin = int(np.where(time == 10.0)[0])
        tmax = int(np.where(time == 18.0)[0])

    else:
        continue

    dstddt = np.zeros(len(std[tmin:tmax]))

    dstddt[0]    = (std[tmin+1] - std[tmin])/dt
    dstddt[1:-1] = (std[(tmin+2):tmax] - std[tmin:(tmax-2)])/(2*dt)
    dstddt[-1]   = (std[tmax] - std[tmax-1])/dt

    avg_dstddt = np.mean(dstddt)

    np.save(path_save+'avg_dstddt_'+runs[i][4:6]+'_'+runs[i][7:]+'.npy',avg_dstddt)

    print('----------------------------------------------------')
    print('For Tperp/Tpar = '+runs[i][5]+', Beta = '+runs[i][8]+':')
    print('dsigma/dt = '+str(avg_dstddt))
    print('From t = '+str(time[tmin])+' s to '+str(time[tmax])+' s')
    print('Total: '+str((tmax-tmin)/2)+' s')
    print(path_save+'avg_dstddt_'+runs[i][4:6]+'_'+runs[i][7:]+'.npy')

print('----------------------------------------------------')

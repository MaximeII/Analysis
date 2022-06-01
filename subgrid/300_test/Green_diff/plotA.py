import numpy as np
import matplotlib.pyplot as plt

path_fig  = '/wrk/users/dubart/analysis/subgrid/300_test/Green_diff/fig/evaluateA/'
path_save = '/wrk/users/dubart/analysis/subgrid/300_test/Green_diff/data/'

dt    = [0.01,0.05,0.001,0.005,0.0001,0.0005]
Dmumu = [0.1,0.05,0.01,0.005,0.001]

dtleg = ['0.01','0.05','0.001','0.005','0.0001','0.0005']
Dleg  = ['0.1','0.05','0.01','0.005','0.001']

fig, axs = plt.subplots(2,1,sharex = True)

axs[0].set_xlim(0.0,12.0)
axs[0].set_ylim(0.29,1.0)
axs[0].set_ylabel('A',fontsize=15)

axs[1].set_xlim(0.0,12.0)
axs[1].set_ylim(0.29,1.0)
axs[1].set_ylabel('A',fontsize=15)
axs[1].set_xlabel('Time (s)',fontsize=15)

for i in range(0,len(dt)):
    A = np.load(path_save+'A_0.29_D0.01_dt'+str(dt[i])+'.npy')    
    dtrange = np.arange(0.0,12.0,dt[i])

    axs[0].plot(dtrange,A,label=dtleg[i])

    print(A)

axs[0].legend(loc='lower right')

for i in range(0,len(Dmumu)):
    A = np.load(path_save+'A_0.29_D'+str(Dmumu[i])+'_dt0.001.npy')
    dtrange = np.arange(0.0,12.0,0.001)    

    axs[1].plot(dtrange,A,label=Dleg[i])

axs[1].legend(loc='lower right')

plt.savefig(path_fig+'convA.png',dpi=300)
print(path_fig+'convA.png')

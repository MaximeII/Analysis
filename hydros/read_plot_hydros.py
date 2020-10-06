import numpy as np
import matplotlib.pyplot as plt
import sys

runs   = ['BCQ','BGA','BCG']
colors = ['tab:blue','tab:red','tab:orange']
dr     = [300,600,900]
d_i    = [132937.78719913246,139580.37722259495,134878.06835706628]

kpar  = np.linspace(0.01,1.5,151)

kpara_AIC = np.zeros([len(runs),len(kpar)])
omega_AIC = np.zeros([len(runs),len(kpar)])
gamma_AIC = np.zeros([len(runs),len(kpar)])

kpara_MM = np.zeros([len(runs),len(kpar)])
omega_MM = np.zeros([len(runs),len(kpar)])
gamma_MM = np.zeros([len(runs),len(kpar)])

fig, axs = plt.subplots(2,1,sharex=True)

for i in range(1,3):
    path_output = '/wrk/users/dubart/analysis/hydros/output/'+runs[i]+'/'
    HYDROS_AIC=np.loadtxt(path_output+runs[i]+'_krange_beta2.log',dtype=float)

    kpara_AIC[i,:] = HYDROS_AIC[:,0]
    omega_AIC[i,:] = HYDROS_AIC[:,2]
    gamma_AIC[i,:] = HYDROS_AIC[:,3]

    axs[0].plot(kpara_AIC[i,:],gamma_AIC[i,:],linewidth=2,color=colors[i],label='$\\Delta r = '+str(dr[i])+'$ km')
    axs[0].axvline(np.pi/(2*dr[i]*1.e3)*d_i[i],0.0,1.0,linewidth=1,linestyle='--',color=colors[i])

for i in range(0,len(runs)-2):
    path_output = '/wrk/users/dubart/analysis/hydros/output/'+runs[i]+'/'
    for j in range(0,len(kpar)):
    
        HYDROS_AIC=np.loadtxt(path_output+runs[i]+'_beta2_'+str(round(kpar[j],2))+'.log',dtype=float)

        kpara_AIC[i,j] = HYDROS_AIC[0]
        omega_AIC[i,j] = HYDROS_AIC[2]
        gamma_AIC[i,j] = HYDROS_AIC[3]

    axs[0].plot(kpara_AIC[i,:],gamma_AIC[i,:],linewidth=2,color=colors[i],label='$\\Delta r = '+str(dr[i])+'$ km')
    axs[0].axvline(np.pi/(2*dr[i]*1.e3)*d_i[i],0.0,1.0,linewidth=1,linestyle='--',color=colors[i])

kpar  = np.linspace(0.01,1.2,121)

for i in range(0,len(runs)):
    path_output = '/wrk/users/dubart/analysis/hydros/output/'+runs[i]+'/'
    for j in range(0,len(kpar)):

        HYDROS_MM=np.loadtxt(path_output+runs[i]+'_beta2_'+str(round(kpar[j],2))+'_45.log',dtype=float)
        
        kpara_MM[i,j] = HYDROS_MM[0]
        omega_MM[i,j] = HYDROS_MM[2]
        gamma_MM[i,j] = HYDROS_MM[3]
 
    axs[1].plot(kpara_MM[i,:],gamma_MM[i,:],linewidth=2,color=colors[i])
    axs[1].axvline(np.pi/(2*dr[i]*1.e3)*d_i[i],0.0,1.0,linewidth=1,linestyle='--',color=colors[i]) 

axs[1].set_xlim(kpar[0],1.5)
axs[0].set_ylim(0.0,0.3)
axs[1].set_ylim(0.0,0.3)
axs[1].set_xlabel('$k * d_i$')
axs[0].legend(loc='best')
axs[0].set_ylabel('$\\gamma / \\omega_{ci}$')
axs[1].set_ylabel('$\\gamma / \\omega_{ci}$')
axs[0].set_title('(a) Proton cyclotron, $\\theta = 0^\circ$')
axs[1].set_title('(b) Mirror, $\\theta = 45^\circ$')

plt.savefig('./AIC_GR_beta2.png',dpi=300)
print('./AIC_GR_beta2.png')

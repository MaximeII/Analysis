import pytools as pt
import numpy as np
from multiprocessing import Pool
import os,sys
from scipy import integrate

numproc = 20 

run = sys.argv[1]
bulkStart  = int(sys.argv[2])
bulkEnd    = int(sys.argv[3])
boundaries = sys.argv[4]

dt = 5.0 # How many seconds to calculate diffusion coeff. Min = 0.5s


if boundaries == 'periodic':

    path_bulk = '/wrk/users/dubart/'+run+'km/bulk/'
    path_save = '/wrk/users/dubart/analysis/subgrid/data/'+run+'/periodic/'
    CellID = int(np.load(path_save+'CellID_box.npy'))

elif boundaries == 'mp':

    path_bulk = '/wrk/users/dubart/'+run+'km/maxwellian_periodic/'
    path_save = '/wrk/users/dubart/analysis/subgrid/data/'+run+'/maxwellian_periodic/'
    CellID = int(np.load(path_save+'CellID_box.npy'))

else:

    path_bulk = '/wrk/group/spacephysics/vlasiator/2D/BCQ/bulk/'
    path_save = '/wrk/users/dubart/analysis/subgrid/data/'+run+'/'
    CellID = int(np.load(path_save+'CellID.npy'))


# Get f(mu,t)
def get_fmu(step):

    if run == '300':
        pt.calculations.pitch_angles_box(vlsvReader=pt.vlsvfile.VlsvReader(path_bulk+'bulk.'+str(step).rjust(7,'0')+'.vlsv'),
                            cellid = CellID, nbins = 32,
                            cosine = True,
                            #plasmaframe = True,
                            outputfile = path_save+'PA/mu'+str(step).rjust(7,'0'))

    elif run == 'BCQ':
        pt.calculations.pitch_angles(vlsvReader=pt.vlsvfile.VlsvReader(path_bulk+'bulk.'+str(step).rjust(7,'0')+'.vlsv'),
                            cellid = CellID, nbins = 32,
                            cosine = True,
                            plasmaframe = True,
                            outputfile = path_save+'PA/mu'+str(step).rjust(7,'0'))

    return

pool = Pool(numproc)
print('Getting fmu ...')
pool.map(get_fmu,range(bulkStart,bulkEnd+1))
pool.terminate()
print('Saved "mu'+str(bulkStart).rjust(7,'0')+'" to "mu'+str(bulkEnd).rjust(7,'0')+'" at '+path_save)


# Get first and second order derivatives in mu
def get_dfdmu_dfdmumu(step):

    lines = []

    filePA = open(path_save+'PA/mu'+str(step).rjust(7,'0'),'r')

    for line in filePA:
        lines.append(line)

    angles   = np.fromstring(lines[3],dtype=float,sep=' ')
    binvalue = np.fromstring(lines[5],dtype=float,sep=' ')

    filePA.close()

    dmu = abs(angles[0]-angles[1])

    dfdmu       = np.zeros(len(binvalue))
    dfdmu[1:-1] = (binvalue[2:] - binvalue[0:-2]) / (2*dmu) # center
    dfdmu[0]    = (binvalue[1]  - binvalue[0])    / dmu     # fw for left edge
    dfdmu[-1]   = (binvalue[-1] - binvalue[-2])   / dmu     # bw for right edge

    dfdmumu       = np.zeros(len(binvalue))
    dfdmumu[1:-1] = (binvalue[2:] - 2.0*binvalue[1:-1] + binvalue[0:-2]) / (dmu*dmu)
    dfdmumu[0]    = (binvalue[2]  - 2.0*binvalue[1]    + binvalue[0])    / (dmu*dmu)
    dfdmumu[-1]   = (binvalue[-1] - 2.0*binvalue[-2]   + binvalue[-3])   / (dmu*dmu)

    return [dfdmu,dfdmumu]

pool = Pool(numproc)
print('Getting dfdmu and dfdmumu ...')
data = pool.map(get_dfdmu_dfdmumu,range(bulkStart,bulkEnd+1))
pool.terminate()
data = np.array(data)

dfdmu   = data[:,0,:]
dfdmumu = data[:,1,:]

np.save(path_save+'dfdmu.npy',dfdmu)
np.save(path_save+'dfdmumu.npy',dfdmumu)

# Get time derivative
def get_dfdt(step):

    if step < (bulkStart+int(dt*2.0)):

        # t = 0
        lines = []

        filePA_prev = open(path_save+'PA/mu'+str(step).rjust(7,'0'),'r')

        for line in filePA_prev:
            lines.append(line)

        angles_prev   = np.fromstring(lines[3],dtype=float,sep=' ')
        binvalue_prev = np.fromstring(lines[5],dtype=float,sep=' ')

        filePA_prev.close()

        # t = 1
        lines = []

        filePA_next = open(path_save+'PA/mu'+str(step+int(dt*2.0)).rjust(7,'0'),'r')

        for line in filePA_next:
            lines.append(line)

        angles_next   = np.fromstring(lines[3],dtype=float,sep=' ')
        binvalue_next = np.fromstring(lines[5],dtype=float,sep=' ')

        filePA_next.close()

        dfdt = (binvalue_next - binvalue_prev) / dt

    elif step > (bulkEnd-int(dt*2.0)):  

        # t = n-1
        lines = []

        filePA_prev = open(path_save+'PA/mu'+str(step-int(dt*2.0)).rjust(7,'0'),'r')

        for line in filePA_prev:
            lines.append(line)

        angles_prev   = np.fromstring(lines[3],dtype=float,sep=' ')
        binvalue_prev = np.fromstring(lines[5],dtype=float,sep=' ')

        filePA_prev.close()

        # t = n
        lines = []

        filePA_next = open(path_save+'PA/mu'+str(step).rjust(7,'0'),'r')

        for line in filePA_next:
            lines.append(line)

        angles_next   = np.fromstring(lines[3],dtype=float,sep=' ')
        binvalue_next = np.fromstring(lines[5],dtype=float,sep=' ')

        filePA_next.close()

        dfdt = (binvalue_next - binvalue_prev) / dt


    else:

        # t-1
        lines = []

        filePA_prev = open(path_save+'PA/mu'+str(step-int(dt*2.0)).rjust(7,'0'),'r')

        for line in filePA_prev:
            lines.append(line)

        angles_prev   = np.fromstring(lines[3],dtype=float,sep=' ')
        binvalue_prev = np.fromstring(lines[5],dtype=float,sep=' ')

        filePA_prev.close()

        # t+1 
        lines = []

        filePA_next = open(path_save+'PA/mu'+str(step+int(dt*2.0)).rjust(7,'0'),'r')

        for line in filePA_next:
            lines.append(line)

        angles_next   = np.fromstring(lines[3],dtype=float,sep=' ')
        binvalue_next = np.fromstring(lines[5],dtype=float,sep=' ')

        filePA_next.close()

        dfdt = (binvalue_next - binvalue_prev)/(2.0*dt)

    return [dfdt]

pool = Pool(numproc)
print('Getting dfdt ...')
data = pool.map(get_dfdt,range(bulkStart,bulkEnd+1))
pool.terminate()
data = np.array(data)

dfdt = data[:,0,:]

np.save(path_save+'dfdt.npy',dfdt)

# Get dmu
lines = []

filePA = open(path_save+'PA/mu'+str(bulkStart).rjust(7,'0'),'r')

for line in filePA:
    lines.append(line)

angles   = np.fromstring(lines[3],dtype=float,sep=' ')

filePA.close()

dmu = abs(angles[0]-angles[1])

# Get diffusion coefficient
def get_Dmumu(step):

    flux = integrate.simps(dfdt[step,:],dx=dmu)     

    Dmumu = - flux / integrate.simps(dfdmumu[step,:],dx=dmu)

    return [Dmumu]

  
pool = Pool(numproc)
print('Getting Dmumu ...')
data = pool.map(get_Dmumu,range(0,dfdt.shape[0]))
pool.terminate()

data = np.array(data)

Dmumu = data[:,0,:]

np.save(path_save+'Dmumu.npy',Dmumu)
print('Saved '+path_save+'Dmumu.npy')

import pytools as pt
import numpy as np
from multiprocessing import Pool
import os,sys


numproc = 20 

run = sys.argv[1]
bulkStart = int(sys.argv[2])
bulkEnd   = int(sys.argv[3])

path_bulk = '/wrk/users/dubart/'+run+'km/Shock/'
#path_bulk = '/wrk/users/dubart/'+run+'km/maxwellian_periodic/'
path_save = '/wrk/users/dubart/'+run+'km/data/'

CellID = int(np.load(path_save+'CellID_Shock.npy'))
#CellID = int(np.load(path_save+'CellID_box.npy'))

# Get f(mu,t)
def get_fmu(step):

    pt.calculations.pitch_angles_box(vlsvReader=pt.vlsvfile.VlsvReader(path_bulk+'bulk.'+str(step).rjust(7,'0')+'.vlsv'),
                        cellid = CellID, nbins = 32,
                        cosine = True,
                        plasmaframe = True,
                        outputfile = path_save+'mu'+str(step).rjust(7,'0'))

    return

pool = Pool(numproc)
print('Getting fmu ...')
pool.map(get_fmu,range(bulkStart,bulkEnd+1))
print('Saved "mu'+str(bulkStart).rjust(7,'0')+'" to "mu'+str(bulkEnd).rjust(7,'0')+'" at '+path_save)


# Get first and second order derivatives in mu
def get_dfdmu_dfdmumu(step):

    lines = []

    filePA = open(path_save+'mu'+str(step).rjust(7,'0'),'r')

    for line in filePA:
        lines.append(line)

    angles   = np.fromstring(lines[3],dtype=float,sep=' ')
    binvalue = np.fromstring(lines[5],dtype=float,sep=' ')

    filePA.close()

    dmu = abs(angles[0]-angles[1])

    dfdmu = (binvalue[2:] - binvalue[0:-2])/(2*dmu)

    dfdmumu = (binvalue[2:] - 2.0*binvalue[1:-1] + binvalue[0:-2])/(dmu*dmu)

    return [dfdmu,dfdmumu]

pool = Pool(numproc)
print('Getting dfdmu and dfdmumu ...')
data = pool.map(get_dfdmu_dfdmumu,range(bulkStart+1,bulkEnd))
data = np.array(data)

dfdmu   = data[:,0,:]
dfdmumu = data[:,1,:]


# Get time derivative
def get_dfdt(step):

    # t-1
    lines = []

    filePA_prev = open(path_save+'mu'+str(step-1).rjust(7,'0'),'r')

    for line in filePA_prev:
        lines.append(line)

    angles_prev   = np.fromstring(lines[3],dtype=float,sep=' ')
    binvalue_prev = np.fromstring(lines[5],dtype=float,sep=' ')

    filePA_prev.close()

    # t+1 
    lines = []

    filePA_next = open(path_save+'mu'+str(step+1).rjust(7,'0'),'r')

    for line in filePA_next:
        lines.append(line)

    angles_next   = np.fromstring(lines[3],dtype=float,sep=' ')
    binvalue_next = np.fromstring(lines[5],dtype=float,sep=' ')

    filePA_next.close()

    dt = 0.5

    dfdt = (binvalue_next[1:-1] - binvalue_prev[1:-1])/(2.0*dt)

    return [dfdt]

pool = Pool(numproc)
print('Getting dfdt ...')
data = pool.map(get_dfdt,range(bulkStart+1,bulkEnd))
data = np.array(data)

dfdt = data[:,0,:]

# Tri Diagonal Matrix Algorithm(a.k.a Thomas algorithm) solver
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


# Get dmu
lines = []

filePA = open(path_save+'mu'+str(bulkStart).rjust(7,'0'),'r')

for line in filePA:
    lines.append(line)

angles   = np.fromstring(lines[3],dtype=float,sep=' ')

filePA.close()

dmu = abs(angles[0]-angles[1])

# Get diffusion coefficient
def get_Dmumu(step):

    Dmumu = TDMAsolver(-dfdmu[step,1:]/(2.0*dmu),dfdmumu[step,:],dfdmu[step,:-1]/(2.0*dmu),dfdt[step,:])

    return [Dmumu]

  
pool = Pool(numproc)
print('Getting Dmumu ...')
data = pool.map(get_Dmumu,range(0,dfdt.shape[0]))
data = np.array(data)

Dmumu = data[:,0,:]


np.save(path_save+'Dmumu_'+run+'_Shock.npy',Dmumu)
print('Saved '+path_save+'Dmumu_'+run+'_Shock.npy')
#np.save(path_save+'Dmumu_'+run+'.npy',Dmumu)
#print('Saved '+path_save+'Dmumu_'+run+'.npy')

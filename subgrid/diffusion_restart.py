import pytools as pt
import numpy as np
import os,sys
from multiprocessing import Pool

numproc = 20 

mode = sys.argv[1] # v or va

dt = 5.0 
t  = 0.5 

path_bulk = '/wrk/group/spacephysics/vlasiator/2D/BCQ/restart/'
path_save = '/wrk/users/dubart/analysis/subgrid/data/BCQ/restart/'

bulkname = 'restart.0001361.vlsv'

RE = 6371e3
x0 = 15.0 * RE
z0 = 10.0 * RE
coord_start = [x0,0.0,z0]

mu0  = 1.26e-6

f = pt.vlsvfile.VlsvReader(path_bulk+bulkname)

CID_start = int(f.get_cellid(coord_start))
if mode == 'v':
    v_start   = f.read_variable('restart_V',CID_start)
elif mode == 'va':
    B        = f.read_variable('b',CID_start)
    Bmag     = np.linalg.norm(B)
    rhom     = f.read_variable('restart_rhom',CID_start)
    v_start  = Bmag / np.sqrt(mu0 *rhom)

print(rhom)

x = x0
v = v_start
z = z0
CellIDs = CID_start

while x > 14.0*RE:

    if mode == 'v':

        x = v[0] * t + x
        z = v[2] * t + z

        CID = int(f.get_cellid([x,0.0,z]))
        v   = f.read_variable('restart_V',CID)

    elif mode == 'va':

        ex = [1.0,0.0,0.0]
 
        theta = np.arccos(np.dot(B,ex) / Bmag)
        
        x = - v*np.cos(theta)*t + x
        z = v*np.sin(theta)*t + z 

        print('vA = '+str(v/1e3)+'km/s')
        print('['+str(x)+',0.0,'+str(z)+']')

        CID = int(f.get_cellid([x,0.0,z]))
        
        B    = f.read_variable('b',CID)
        Bmag = np.linalg.norm(B)
        rhom = f.read_variable('restart_rhom',CID)
        v    = Bmag / np.sqrt(mu0 *rhom) 

    CellIDs = np.append(CellIDs,CID)

print('Got CellIDs for '+mode)
np.save(path_save+'CellIDs_'+mode+'.npy',CellIDs)

# Get f(mu,t)
def get_fmu(step):
#for i in range(0,len(CellIDs)):

    pt.calculations.pitch_angles(vlsvReader=pt.vlsvfile.VlsvReader(path_bulk+bulkname),
                        cellid = CellIDs[step], nbins = 32,
                        cosine = True,
                        plasmaframe = True,
                        outputfile = path_save+'PA/mu'+str(step).rjust(4,'0')+'_'+mode)

#    print('Saved mu'+str(i).rjust(4,'0')+'_'+mode)
    return

pool = Pool(numproc)
print('Getting f(mu) ...')
pool.map(get_fmu,range(0,len(CellIDs)))
pool.terminate()
print('Saved f(mu) at '+path_save+'PA/')

# Get first and second order derivatives in mu
def get_dfdmu_dfdmumu(step):

    lines = []

    filePA = open(path_save+'PA/mu'+str(step).rjust(4,'0')+'_'+mode,'r')

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
data = pool.map(get_dfdmu_dfdmumu,range(0,len(CellIDs)))
pool.terminate()
data = np.array(data)

dfdmu   = data[:,0,:]
dfdmumu = data[:,1,:]

np.save(path_save+'dfdmu_'+mode+'.npy',dfdmu)
np.save(path_save+'dfdmumu_'+mode+'.npy',dfdmumu)

# Get time derivative
def get_dfdt(step):

    if step < (0 + int(dt*2.0)):

        # t = 0
        lines = []

        filePA_prev = open(path_save+'PA/mu'+str(step).rjust(4,'0')+'_'+mode,'r')

        for line in filePA_prev:
            lines.append(line)

        angles_prev   = np.fromstring(lines[3],dtype=float,sep=' ')
        binvalue_prev = np.fromstring(lines[5],dtype=float,sep=' ')

        filePA_prev.close()

        # t = 1
        lines = []

        filePA_next = open(path_save+'PA/mu'+str(step+int(dt*2.0)).rjust(4,'0')+'_'+mode,'r')

        for line in filePA_next:
            lines.append(line)

        angles_next   = np.fromstring(lines[3],dtype=float,sep=' ')
        binvalue_next = np.fromstring(lines[5],dtype=float,sep=' ')

        filePA_next.close()

        dfdt = (binvalue_next - binvalue_prev) / dt

    elif step > (len(CellIDs)-1 - int(dt*2.0)):

        # t = n-1
        lines = []

        filePA_prev = open(path_save+'PA/mu'+str(step-int(dt*2.0)).rjust(4,'0')+'_'+mode,'r')

        for line in filePA_prev:
            lines.append(line)

        angles_prev   = np.fromstring(lines[3],dtype=float,sep=' ')
        binvalue_prev = np.fromstring(lines[5],dtype=float,sep=' ')

        filePA_prev.close()

        # t = n
        lines = []

        filePA_next = open(path_save+'PA/mu'+str(step).rjust(4,'0')+'_'+mode,'r')

        for line in filePA_next:
            lines.append(line)

        angles_next   = np.fromstring(lines[3],dtype=float,sep=' ')
        binvalue_next = np.fromstring(lines[5],dtype=float,sep=' ')

        filePA_next.close()

        dfdt = (binvalue_next - binvalue_prev) / dt


    else:

        # t-1
        lines = []

        filePA_prev = open(path_save+'PA/mu'+str(step-int(dt*2.0)).rjust(4,'0')+'_'+mode,'r')

        for line in filePA_prev:
            lines.append(line)

        angles_prev   = np.fromstring(lines[3],dtype=float,sep=' ')
        binvalue_prev = np.fromstring(lines[5],dtype=float,sep=' ')

        filePA_prev.close()

        # t+1 
        lines = []

        filePA_next = open(path_save+'PA/mu'+str(step+int(dt*2.0)).rjust(4,'0')+'_'+mode,'r')

        for line in filePA_next:
            lines.append(line)

        angles_next   = np.fromstring(lines[3],dtype=float,sep=' ')
        binvalue_next = np.fromstring(lines[5],dtype=float,sep=' ')

        filePA_next.close()

        dfdt = (binvalue_next - binvalue_prev)/(2.0*dt)

    return [dfdt]

pool = Pool(numproc)
print('Getting dfdt ...')
data = pool.map(get_dfdt,range(0,len(CellIDs)))
pool.terminate()
data = np.array(data)

dfdt = data[:,0,:]

np.save(path_save+'dfdt_'+mode+'.npy',dfdt)

# Get dmu
lines = []

filePA = open(path_save+'PA/mu0000_'+mode,'r')

for line in filePA:
    lines.append(line)

angles   = np.fromstring(lines[3],dtype=float,sep=' ')

filePA.close()

dmu = abs(angles[0]-angles[1])

# Get diffusion coefficient
def get_Dmumu(step):

    Mat = np.zeros([dfdmumu.shape[1],dfdmumu.shape[1]])

    for j in range(0,dfdmumu.shape[1]):
        for k in range(0,dfdmumu.shape[1]):

            if k == j:
                Mat[j,k] = dfdmumu[step,j] # Diagonal
            elif k == j+1:
                Mat[j,k] = dfdmu[step,j]/(2*dmu) # upper
            elif k == j-1:
                Mat[j,k] = -dfdmu[step,j]/(2*dmu) # lower


#    Mat[0,-1] = -dfdmu[step,0]/(2*dmu) # added for periodic boundary conditions
#    Mat[-1,0] = dfdmu[step,-1]/(2*dmu)

    Dmumu = np.linalg.solve(Mat,dfdt[step,:])

    return [Dmumu]


pool = Pool(numproc)
print('Getting Dmumu ...')
data = pool.map(get_Dmumu,range(0,dfdt.shape[0]))
pool.terminate()

data = np.array(data)

Dmumu = data[:,0,:]

np.save(path_save+'Dmumu_'+mode+'.npy',Dmumu)
print('Saved '+path_save+'Dmumu_'+mode+'.npy')


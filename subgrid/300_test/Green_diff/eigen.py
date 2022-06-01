import numpy as np
import sys

run = sys.argv[1]

path_bulk = '/wrk/users/dubart/diff_test/dmumu/'+run+'/bulk/'
path_save = '/wrk/users/dubart/analysis/subgrid/300_test/Green_diff/data/mu/'+run+'/'
path_fig  = '/wrk/users/dubart/analysis/subgrid/300_test/Green_diff/fig/mu/'+run+'/'

Matrix = np.load(path_save+'Matrix_'+run+'.npy')

for i in range(0,Matrix.shape[2]):
    eigvalues, eigvector = np.linalg.eig(Matrix[:,:,i])
    
    print('------ i = '+str(i)+' ------')
    print('If no zero eigenvalues, then TRUE: '+str(np.all(eigvalues)))
    print(Matrix[:,:,i])
    #print('EIGENVALUES')
    #print(eigvalues)
    #print('EIGENVECTOR') 
    #print(eigvector)

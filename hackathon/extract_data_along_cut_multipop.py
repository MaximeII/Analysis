import pytools as pt
import numpy as np
import sys
from variable import get_data, get_name, get_units
from multiprocessing import Pool

numproc = 40

path_bulk = '/wrk/group/spacephysics/vlasiator/3D/EGI/bulk/dense_cold_hall1e5_afterRestart374/'
path_save = '/wrk/users/dubart/analysis/hackathon/data/'

bulkStart = 662
bulkEnd   = 1503

RE = 6371e3

xmin = -25.0 *RE
xmax = -8.0  *RE
zmin = 0.0   *RE
zmax = 0.0   *RE
ymin = -4.0   *RE
ymax = -4.0   *RE
#xmin = -7.04e8 + 16e6
#ymin = -3.68e8 + 16e6
#zmin = -3.68e8 + 16e6


point1 = [xmin,ymin,zmin]
point2 = [xmax,ymax,zmax]

pop = 'proton'

filetitle = 'EGI_cut_x_y0_z0'

def cuts(step):

    filename = 'bulk1.'+str(step).rjust(7,'0')+'.vlsv'

    #open a vlsv file

    print(filename)

    f = pt.vlsvfile.VlsvReader(path_bulk+filename)

    #cut = pt.calculations.cut_through(f,point1,point2)

    #print(f.get_cellid(point1))
    #print(f.get_cellid(point2))

    cellID1 = int(f.get_cellid(point1))
    cellID2 = int(f.get_cellid(point2))

    cut = np.arange(cellID1,cellID2,1)

    #prof_re = get_data(cut[1])/RE
    
    prof_re = np.zeros(len(cut))
    for i in range(0,len(cut)):
        prof_re[i] = f.get_cell_coordinates(cut[i])[0]

    prof_start = prof_re[0]
    prof_stop  = prof_re[-1]

    pr_np    = get_data(f.read_variable(pop+"/vg_rho",cut))
    pr_Vpx   = get_data(f.read_variable_info(pop+"/vg_v",cut,"x"))
    pr_Vpy   = get_data(f.read_variable_info(pop+"/vg_v",cut,"y"))
    pr_Vpz   = get_data(f.read_variable_info(pop+"/vg_v",cut,"z"))
    pr_Tpara = get_data(f.read_variable_info(pop+"/vg_t_parallel",cut))
    pr_Tperp = get_data(f.read_variable_info(pop+"/vg_t_perpendicular",cut))

    pr_Bmag = get_data(f.read_variable_info("vg_b_vol",cut,"magnitude"))
    pr_Bx   = get_data(f.read_variable_info("vg_b_vol",cut,"x"))
    pr_By   = get_data(f.read_variable_info("vg_b_vol",cut,"y"))
    pr_Bz   = get_data(f.read_variable_info("vg_b_vol",cut,"z"))

    #pr_Emag = get_data(f.read_variable_info("vg_e_vol",cut[0].data,"magnitude"))
    #pr_Ex   = get_data(f.read_variable_info("vg_e_vol",cut[0].data,"x"))
    #pr_Ey   = get_data(f.read_variable_info("vg_e_vol",cut[0].data,"y"))
    #pr_Ez   = get_data(f.read_variable_info("vg_e_vol",cut[0].data,"z"))

    print("Writing data...")

    f_cut = open(path_save+filetitle+'_'+str(step).rjust(7,'0')+'.dat',"w")

    f_cut.write("x y z Bmag Bx By Bz np Vpx Vpy Vpz Tpara Tperp\n")

    for i in range(0,len(prof_re)):
        f_cut.write(str(prof_re[i])+" "+str(ymin)+" "+str(zmin)+" "+str(pr_Bmag[i])+" "+str(pr_Bx[i])+" "+str(pr_By[i])+" "+str(pr_Bz[i])+" "+str(pr_np[i])+" "+str(pr_Vpx[i])+" "+str(pr_Vpy[i])+" "+str(pr_Vpz[i])+" "+str(pr_Tpara[i])+" "+str(pr_Tperp[i])+"\n")

    f_cut.close()

    return

pool = Pool(numproc)
pool.map(cuts,range(bulkStart,bulkEnd+1))
pool.terminate()


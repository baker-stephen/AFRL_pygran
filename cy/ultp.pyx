import numpy as np
import matplotlib.pyplot as plt
import time

import sys
sys.path.insert(1,'../')
from param_defn import PD

def go(params: PD, name_mod=""):

    print("Inputs received:")
    wrt_dir = params.out_dir
    print(wrt_dir)
    Dp = params.DP
    print(Dp)
    ID = params.ID
    print(ID)
    z_0 = params.z0()
    print(z_0)
    z_M = params.zM()
    print(z_M)
    z0_cutoff = params.z0_cutoff()
    print(z0_cutoff)
    zM_cutoff = params.zM_cutoff()
    print(zM_cutoff)
    res = params.res()
    print(res)

    start_time = time.time()
    loaded = []

    cdef double pipe_in_rad = ID / 2.0

    cdef int x_res = res[0]
    cdef double x0 = -pipe_in_rad
    cdef double xM = pipe_in_rad
    cdef double dx = (xM - x0) / (x_res + 1.0)

    cdef int y_res = res[1]
    cdef double y0 = -pipe_in_rad
    cdef double yM = pipe_in_rad
    cdef double dy = (yM - y0) / (y_res + 1.0)

    cdef int z_res = res[2]
    cdef double z0 = z_0*1.0
    cdef double zM = z_M*1.0
    cdef double dz = (zM - z0) / (z_res + 1.0)

    # Load the previously generated 3D grid
    with open(wrt_dir+'filled_array_'+name_mod+'_'+str(x_res)+'-'+str(y_res)+'-'+str(z_res)+'.npy', 'rb') as f:
        loaded = np.load(f)
        f.close()

    #Calculate average porosity as a function of z position
    avg_z = []
    cdef double dV = dx*dy*dz #Volume element
    cdef double total_vol_z = np.pi * pipe_in_rad ** 2.0 * dz #Total volume of "disk" with height dz
    cdef double sumz,fill_vol_z,p_z = 0.0 #Cython initialize variables for small speedup
    for z in range(0,z_res):
        sumz = 0.0
        # At this z postion, sum all contributions in the x and y directions
        for x in range(x_res+1):
            for y in range(y_res+1):
                sumz += loaded[x][y][z]
        if z%100==0:
            print("completed z =",z)
        fill_vol_z = sumz * dV #Gives the physical volume for this disk
        p_z = (total_vol_z-fill_vol_z)/total_vol_z #Calculate porosity
        avg_z.append(p_z)

    print("finished z, writing out avg_z numpy array")
    with open(wrt_dir+'avg_z_'+name_mod+'_'+str(x_res)+'-'+str(y_res)+'-'+str(z_res)+'.npy', 'wb') as f:
        np.save(f, np.array(avg_z))
    print("finished writing np z array")
    print("mean poros for zs:",np.mean(avg_z))
    print("std dev in z:",np.std(avg_z))

    figz, axz = plt.subplots()
    axz.plot(np.linspace(z0,zM,len(avg_z)),avg_z)
    plt.savefig(wrt_dir+'z-poros-data_'+name_mod+str(x_res)+'-'+str(y_res)+'-'+str(z_res)+'.png')

    print("--- finish z in %s seconds ---" % (time.time()-start_time))

    #Begin average porosity as a function of radius

    #Initialize array to write to
    cdef int r_res = res[3]
    avg_rs = np.zeros(r_res+1) #Porosity values will be stored here
    rs = np.linspace(0,pipe_in_rad,len(avg_rs)) #This is the 1D radius "grid"
    cdef double dr = pipe_in_rad/(r_res+1.0)
    cdef double pr_sq = pipe_in_rad**2.0

    zstart = int(round((z0_cutoff-z0)/dz)) #Logical index of where we start calculating porosity
    print("zstart:",zstart)
    zend = int(round((zM_cutoff-z0)/dz)) #Logical index of where we stop calculating porosity
    print("znd:",zend)
    for z in range(zstart,zend):

        #Every 10% print progress
        if (z-zstart)%int((zend-zstart)*0.1)==0:
            print(str(round(100*(z-zstart)/(zend-zstart))) + "%" +" complete")

        for x in range(x_res+1):
            x_phys = abs(x*dx+x0) #physical coordinate of the center of this node
            if abs(x_phys) >= pipe_in_rad+dx:
                #This x is not within the domain
                continue
            elif abs(x_phys) >= pipe_in_rad:
                #This x is on the edge of the domain
                avg_rs[r_res] += loaded[x][int(y_res/2)][z]
                continue
            #Calculate the possible range of y values for this x
            y_min = -np.sqrt(pr_sq-x_phys**2)
            y_max = -y_min
            y_min_li = max(round(-.5 + ((y_min - y0) / dy)), 0)
            y_max_li = min(round(.5 + ((y_max - y0) / dy)), y_res)

            #Iterate over possible y indices
            for y in range(y_min_li,y_max_li):
                y_phys = abs(y*dy+y0)
                if abs(y_phys) >= pipe_in_rad+dy:
                    continue
                elif abs(y_phys) >= pipe_in_rad:
                    avg_rs[r_res] += loaded[x][y][z]
                    continue
                r_phys = np.sqrt(x_phys**2+y_phys**2)
                #now scatter to radii nodes outward and inward (linear extrapolation)
                r_node_in = 0
                r_node_out = r_res
                if r_phys==0.0:
                    avg_rs[0] += loaded[x][y][z]
                elif r_phys>=pipe_in_rad+dr:
                    continue
                elif r_phys>=pipe_in_rad:
                    avg_rs[r_res] += loaded[x][y][z]
                else:
                    r_node_in = int(r_phys/dr)
                    if r_node_in >= r_res:
                        #Prevent out-of-bounds
                        continue

                    r_node_out = r_node_in+1

                    delta = r_phys/dr-r_node_in
                    avg_rs[r_node_in] += (1-delta)*loaded[x][y][z]
                    avg_rs[r_node_out] += delta * loaded[x][y][z]


    #Compute annular volume element for each node to determine void volume fraction:
    for r in range(0,r_res+1):
        r_in = (r - .5) * dr
        if r == 0:
            r_in = 0
        r_out = (r+.5)*dr
        cyl_out = np.pi*r_out**2*(zM_cutoff-z0_cutoff)
        cyl_in = np.pi * r_in ** 2 * (zM_cutoff - z0_cutoff)
        total_vol = cyl_out-cyl_in
        dV = dx*dy*dz
        fill_vol_r = avg_rs[r]*dV
        avg_rs[r] = (total_vol-fill_vol_r)/total_vol
    print("finished everything, writing out avg_rs numpy array")
    with open(wrt_dir+'avg_r_'+name_mod+'_'+str(x_res)+'-'+str(y_res)+'-'+str(z_res)+'-'+str(r_res)+'.npy', 'wb') as f:
        np.save(f, np.array(avg_rs))
    print("finished writing np array")
    print("mean poros for rs:",np.mean(avg_rs))
    print("std dev in r:",np.std(avg_rs))
    figr, axr = plt.subplots()
    axr.plot(rs,avg_rs)
    plt.savefig(wrt_dir+'r-poros-data_'+name_mod+str(x_res)+'-'+str(y_res)+'-'+str(z_res)+'-'+str(r_res)+'.png')
    rt = time.time() - start_time
    print("--- %s seconds ---" % rt)



    with open(wrt_dir+'outputs_'+name_mod+'.txt', 'a') as out:
        out.write("mean poros for zs:" + str(np.mean(avg_z)) + '\n')
        out.write("std dev in z:" + str(np.std(avg_z)) + '\n')
        out.write("mean poros for rs:" + str(np.mean(avg_rs)) +'\n')
        out.write("std dev in r:" + str(np.std(avg_rs)) + '\n')
        out.write('ultimate_poros runtime: ' + str(rt) + '\n')
        out.close()

    print("Complete radial porosity determination")

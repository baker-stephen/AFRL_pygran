import numpy as np
import time

import sys
sys.path.insert(1,'../')
from param_defn import PD

#TODO: document

def go(params: PD):

    print("Beginning overall porosity determination and grid construction.")
    print("Parameters received:")
    final_step_out = params.final_step() - (params.final_step() % 50000)
    import_csv = params.out_dir+'csvs/particles'+str(final_step_out)+'.csv'
    print(import_csv)
    Dp = params.DP
    print(Dp)
    ID = params.ID
    print(ID)
    z_0 = params.z0()
    print(z_0)
    z_M = params.zM()
    print(z_M)
    res = params.res()
    print(res)

    print("start")
    start_time = time.time()

    positions = []

    with open(import_csv, 'r') as r:
        r.readline()
        for line in r:
            xyz = [float(p) for p in line.split(",")]
            if xyz[2]<0:
                #This is a placeholder particle for Paraview, not part of the simulation. Exclude it.
                continue
            positions.append(xyz)
        r.close()

    print("opened", import_csv)

    wrt_dir = import_csv[:import_csv.rindex("/")-4]

    print("writing to directory: "+ wrt_dir)

    positions = sorted(positions, key=lambda x: x[2])

    cdef double ball_r = Dp/2.0
    print("ball_r:",ball_r)
    cdef double ball_r_sq = ball_r**2.0 #precompute r^2
    print("ball_rsq:", ball_r_sq)
    cdef double full_ball_vol = (4.0 / 3.0) * np.pi * ball_r ** 3.0
    cdef double total_balls_vol = 0.0

    cdef double pipe_in_rad = ID / 2.0

    x_res = res[0]
    cdef double x0 = -pipe_in_rad
    cdef double xM = pipe_in_rad
    cdef double dx = (xM - x0) / (x_res + 1.0)

    y_res = res[1]
    cdef double y0 = -pipe_in_rad
    cdef double yM = pipe_in_rad
    cdef double dy = (yM - y0) / (y_res + 1.0)

    z_res = res[2]
    cdef double z0 = z_0*1.0 #TODO adjust z0,zM
    cdef double zM = z_M*1.0
    print('Height in consideration: %d' % (zM-z0))
    cdef double dz = (zM - z0) / (z_res + 1.0)

    is_filled = []
    for i in range(x_res + 1):
        is_filled.append([])
        for j in range(y_res + 1):
            is_filled[i].append(np.zeros(z_res + 1))
    print("begin main loop")
    cdef double xc = 0.0
    cdef double yc = 0.0
    cdef double zc = 0.0
    cdef double h, vol_cap = 0.0
    for p_num,p in enumerate(positions):
        if p_num%int(len(positions)*0.1)==0:
            print("%d".format(round(100*p_num/len(positions)))+'%'+ ' complete')
        xc = p[0]
        yc = p[1]
        zc = p[2]

        # part of the sphere is outside the standard pipe radius, no longer valid
        if np.sqrt(xc**2+yc**2) + ball_r - 0.08 > pipe_in_rad: #TODO: -0.08 needed for 1.59", check for 1.029" (not necessary 0.26")
            print("invalid sphere at z:",p[2])
            print("radius:",np.sqrt(xc**2+yc**2) + ball_r)
            print('pipe rad:',pipe_in_rad)
            print('diff:',np.sqrt(xc**2+yc**2)-pipe_in_rad)
            raise Exception("zM too high. Keep below " + str(zc-ball_r))

        # If the entirety of the current position is greater than height max, we're done, since sorted by z
        if zc - ball_r > zM:
            break

        # the entirety of the current ball is below the height minimum, disregard
        if zc + ball_r < z0:
            continue

        # The whole ball is within the region of interest
        elif zc + ball_r < zM and zc - ball_r > z0:
            total_balls_vol += full_ball_vol

        # A portion of the ball is within the region of interest, subtract the volume that is above
        elif zc - ball_r > z0:
            h = zc + ball_r - zM
            vol_cap = (1.0 / 3.0) * np.pi * h ** 2.0 * (3.0 * ball_r - h)
            total_balls_vol += full_ball_vol - vol_cap

        # A portion of the ball is within the region of interest, subtract the volume that is below
        else:
            h = zc + ball_r - z0
            vol_cap = (1.0 / 3.0) * np.pi * h ** 2.0 * (3.0 * ball_r - h)
            total_balls_vol += vol_cap

        #Calculate porosity as a function of x,y,z

        # lowest extent of this ball
        z_low = max(z0, zc - ball_r)
        # nearest logical index to the lowest extent of this ball
        z_low_li = max(round(-.5 + ((z_low - z0) / dz)), 0)
        # repeat for highest extent
        z_high = min(zM, zc + ball_r)
        z_high_li = min(round(.5 + ((z_high - z0) / dz)), z_res)
        z_limit = True #True when z distance from center of ball = ball radius
        z_frac = 1.0 #How much of the ball covers this z logical coordinate
        for z in range(z_low_li, z_high_li):
            # given this z coordinate, calculate the x extents:
            z_phys = (z * dz + z0)  # recover physical z
            z_ball = abs(z_phys - zc)  # logical z distance from center of ball
            z_ball_sq = z_ball ** 2
            x_left_li = 0
            x_right_li = 0
            x_left = xc
            x_right = xc
            xr = 0
            if z_ball > ball_r:
                z_frac = 1-((z_ball - ball_r) / dz)
                if z_frac > 1.0:
                    print("z_frac bigger:",z_frac)
                    break
            else:
                z_frac = 1.0
            if z_ball_sq >= ball_r_sq:
                z_limit = True
                if z_frac > 1.0:
                    print("z_frac bigger:",z_frac)
                    break
                x_left_li = round((xc - x0) / dx)
                x_right_li = x_left_li+1
            else:
                xr = np.sqrt(ball_r_sq - z_ball_sq)
                z_limit = False
                x_left = max(x0, xc - xr)
                x_left_li = max(round(-.5+((x_left - x0) / dx)), 0)
                x_right = min(xM, xc + xr)
                x_right_li = min(round(.5+((x_right - x0) / dx)), x_res)
            x_frac = 1.0
            for x in range(x_left_li, x_right_li):
                x_phys = (x * dx + x0)
                x_ball = abs((x_phys - x_left) - xr) #changed from ball_r to xr
                x_ball_sq = x_ball**2
                y_back = yc
                y_back_li = 0
                y_front = yc
                y_front_li = 0
                yr = 0
                if x_ball > xr:
                    x_frac = (x_ball - xr) / dx
                    if x_frac > 1.0:
                        print("x_frac bigger:", x_frac)
                        break
                else:
                    x_frac = 1.0

                if z_limit or z_ball_sq + x_ball_sq >= ball_r_sq:
                    y_back_li = round((yc - y0) / dy)
                    y_front_li = y_back_li+1
                else:
                    yr = np.sqrt(ball_r_sq - z_ball_sq - x_ball_sq)
                    y_back = max(y0, yc - yr)
                    y_back_li = max(round(-.5+((y_back - y0) / dy)), 0)
                    y_front = min(yM, yc + yr)
                    y_front_li = min(round(.5+((y_front - y0) / dy)), y_res)
                for y in range(y_back_li, y_front_li):
                    y_phys = (y * dy + y0)
                    y_ball = abs((y_phys - y_back) - yr) #changed
                    y_frac = 1.0
                    if y_ball > yr:
                        y_frac = 1-((y_ball - yr) / dy)
                        if y_frac < 0.0:
                            print("y_frac negative:", y_frac)
                            y_frac = 0.0
                            break
                    else:
                        y_frac = 1.0
                    if z_frac*x_frac*y_frac>1.0:
                        print("filling with too big:",z_frac*x_frac*y_frac)
                    elif z_frac*x_frac*y_frac<0:
                        print("filling with negative:", z_frac * x_frac * y_frac)
                    is_filled[x][y][z] = (z_frac*x_frac*y_frac)**(1/3)

    print("finished filling, writing out numpy array")
    #TODO: adjust output array name
    with open(wrt_dir+'filled_array_'+str(x_res)+'-'+str(y_res)+'-'+str(z_res)+'.npy', 'wb') as f:
        np.save(f, np.array(is_filled))
    print("finished writing np array")

    pipe_vol = np.pi * pipe_in_rad ** 2 * (zM - z0)

    porosity = (pipe_vol - total_balls_vol) / pipe_vol

    print("total porosity:", porosity)

    sum_zs = []
    sum_y_zs = []
    for i in range(x_res+1):
        sum_zs.append([])
        for j in range(y_res+1):
            sum_zs[i].append(np.sum(is_filled[i][j]))
        sum_y_zs.append(np.sum(sum_zs[i]))
    calc_sum = np.sum(sum_y_zs)
    fill_vol = calc_sum * dx * dy * dz

    print("pipe vol:", pipe_vol)
    summed_porous = (pipe_vol - fill_vol) / pipe_vol
    print("calc from sum:", summed_porous)

    err = abs(summed_porous-porosity)/summed_porous
    print("err:", err)

    rt = time.time() - start_time
    print("--- %s seconds ---" % rt)

    print("writing outputs.txt")

    with open(wrt_dir+'outputs.txt', 'a') as out:
        out.write('\ntotal porosity: ' + str(porosity) + '\n\n')
        out.write(str(x_res)+'-'+str(y_res)+'-'+str(z_res)+'\n')
        out.write("z0="+str(z_0)+",zM="+str(z_M)+'\n')
        out.write('calc from sum: ' + str(summed_porous) + '\n')
        out.write('err: ' + str(err) + '\n')
        out.write('porosity_determination_new runtime: ' + str(rt) + '\n')
        out.close()

    is_filled.clear()

    print("Done!")




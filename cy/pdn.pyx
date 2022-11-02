import numpy as np
import time

import sys
sys.path.insert(1,'../')
from param_defn import PD

def go(params: PD):

    print("Beginning overall porosity determination and grid construction.")
    print("Parameters received:")
    #Since output is only written every 50000 steps, the final recorded position of the particles can be found via:
    final_step_out = params.final_step() - (params.final_step() % 50000)
    import_csv = '../'+params.out_dir+'particles'+str(final_step_out)+'.csv'
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

    #Read the csv array and record the x,y,z position of each particle

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

    # Write the rest of the files to the parent folder (above /csvs/)
    wrt_dir = import_csv[:import_csv.rindex("/")-4]

    print("writing to directory: "+ wrt_dir)

    #Sort positions by z
    positions = sorted(positions, key=lambda x: x[2])

    # Precompute some values to reduce computation time
    cdef double ball_r = Dp/2.0
    print("ball_r:",ball_r)
    cdef double ball_r_sq = ball_r**2.0
    print("ball_rsq:", ball_r_sq)
    cdef double full_ball_vol = (4.0 / 3.0) * np.pi * ball_r ** 3.0
    cdef double total_balls_vol = 0.0

    cdef double pipe_in_rad = ID / 2.0

    #Define 3D grid dimensions

    x_res = res[0]
    cdef double x0 = -pipe_in_rad
    cdef double xM = pipe_in_rad
    cdef double dx = (xM - x0) / (x_res + 1.0)

    y_res = res[1]
    cdef double y0 = -pipe_in_rad
    cdef double yM = pipe_in_rad
    cdef double dy = (yM - y0) / (y_res + 1.0)

    z_res = res[2]
    cdef double z0 = z_0*1.0 #Make sure Cython interprets these values as doubles
    cdef double zM = z_M*1.0
    print('Height in consideration: %d' % (zM-z0))
    cdef double dz = (zM - z0) / (z_res + 1.0)

    # is_filled is a three-dimensional array describing whether or not a sphere overlaps that point.
    # i.e. is_filled[x][y][z] = 1 if (x,y,z) is within a sphere, 0 if not, and if the edge of the sphere falls within
    # half a step (dx/dy/dz) of (x,y,z), it will be a fraction representing the relative proximity of that edge to (x,y,z).
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
    # Iterate through all sphere positions
    for p_num,p in enumerate(positions):
        # Progress printing
        if p_num%int(len(positions)*0.1)==0:
            print(str(round(100*p_num/len(positions)))+'%'+ ' complete')

        # These are the center of sphere p
        xc = p[0]
        yc = p[1]
        zc = p[2]

        # If the entirety of the current position is greater than height max, we're done, since sorted by z
        if zc - ball_r > zM:
            break

        # the entirety of the current ball is below the height minimum, disregard
        if zc + ball_r < z0:
            continue

        # Determine whether part of the sphere is outside the standard pipe radius, and is therefore not valid.
        # This can happen in a multitude of ways. The current sphere may be sitting in the funnel section of the
        # simulation. Or, if the funnel used to fill the sphere employed rounded edges, the rounding may have increased
        # the inner diameter of the pipe within the specified region of interest. To fix either of these, change the
        # second value in the post_zs dictionary (in the param_defn.py file) for your inner diameter, decreasing it
        # below what is specified in the exception.
        if np.sqrt(xc ** 2 + yc ** 2) + ball_r - 0.08 > pipe_in_rad:
            #TODO: I have included the extra -0.08 buffer which is needed for 1.59" ID, but it is not necessary 0.26".
            #   This buffer is needed since the stl file is not a perfectly round cylinder, and has some small edges
            #   which may protrude out slightly but are otherwise perfectly valid. Steps need to be taken to reduce
            #   this value in every stl file. Oversimplification in Meshlab may be to blame.
            print("invalid sphere at z:", p[2])
            print("radius:", np.sqrt(xc ** 2 + yc ** 2) + ball_r)
            print('pipe rad:', pipe_in_rad)
            print('diff:', np.sqrt(xc ** 2 + yc ** 2) - pipe_in_rad)
            raise Exception("zM too high. Keep below " + str(zc - ball_r))

        # First, calculate porosity by simply adding the volume of each sphere within the region of interest.

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

        # Now, calculate porosity as a function of x,y,z
        # To do this, we iterate between the minimum and maximum z extents of the sphere. For each value of z, we
        # determine the maximum x and y extents, and iterate through those to fill the is_filled array with all the
        # (x,y,z) positions that this sphere occupies.

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
            x_left_li = 0 #leftmost logical x index
            x_right_li = 0
            x_left = xc
            x_right = xc
            xr = 0
            if z_ball > ball_r:
                z_frac = 1-((z_ball - ball_r) / dz)
                if z_frac > 1.0:
                    #There's been some sort of mistake, this value should never > 1
                    print("z_frac too large:",z_frac)
                    break
            else:
                z_frac = 1.0
            if z_ball_sq >= ball_r_sq:
                # The z position is at the edge of the sphere. Don't need to scan in the x and y directions anymore.
                z_limit = True
                if z_frac > 1.0:
                    print("z_frac too large:",z_frac)
                    break
                x_left_li = round((xc - x0) / dx)
                x_right_li = x_left_li+1
            else:
                # xr is the radius of the circle that is the intersection of the x-y plane at this z value and sphere p.
                xr = np.sqrt(ball_r_sq - z_ball_sq)
                z_limit = False
                x_left = max(x0, xc - xr) # use max() so we don't scan outside of simulation bounds.
                x_left_li = max(round(-.5+((x_left - x0) / dx)), 0)
                x_right = min(xM, xc + xr)
                x_right_li = min(round(.5+((x_right - x0) / dx)), x_res)
            x_frac = 1.0
            for x in range(x_left_li, x_right_li):
                # Repeat the process for x
                x_phys = (x * dx + x0)
                x_ball = abs((x_phys - x_left) - xr)
                x_ball_sq = x_ball**2
                y_back = yc
                y_back_li = 0
                y_front = yc
                y_front_li = 0
                yr = 0
                if x_ball > xr:
                    x_frac = (x_ball - xr) / dx
                    if x_frac > 1.0:
                        print("x_frac too large:", x_frac)
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
                    # Repeat again for y
                    y_phys = (y * dy + y0)
                    y_ball = abs((y_phys - y_back) - yr)
                    y_frac = 1.0
                    if y_ball > yr:
                        y_frac = 1-((y_ball - yr) / dy)
                        if y_frac < 0.0:
                            # Remnant from debugging
                            print("y_frac negative:", y_frac)
                            y_frac = 0.0
                            break
                    else:
                        y_frac = 1.0
                    if z_frac*x_frac*y_frac>1.0:
                        print("filling with too big:",z_frac*x_frac*y_frac)
                    elif z_frac*x_frac*y_frac<0:
                        print("filling with negative:", z_frac * x_frac * y_frac)
                    # Finally, we can write this x,y,z position
                    is_filled[x][y][z] = (z_frac*x_frac*y_frac)**(1/3)

    # Save this as a binary numpy array file for use in the next step, and easy recovery later.
    print("finished filling, writing out numpy array")
    with open(wrt_dir+'filled_array_'+str(x_res)+'-'+str(y_res)+'-'+str(z_res)+'.npy', 'wb') as f:
        np.save(f, np.array(is_filled))
    print("finished writing np array")

    #Final porosity calculations:

    pipe_vol = np.pi * pipe_in_rad ** 2 * (zM - z0) # Total volume to be filled

    print("pipe vol:", pipe_vol)

    porosity = (pipe_vol - total_balls_vol) / pipe_vol # Simplest, most accurate porosity determination method.

    print("total porosity:", porosity)

    # Calculate porosity by summing the is_filled array
    sum_zs = []
    sum_y_zs = []
    for i in range(x_res+1):
        sum_zs.append([])
        for j in range(y_res+1):
            sum_zs[i].append(np.sum(is_filled[i][j]))
        sum_y_zs.append(np.sum(sum_zs[i]))
    calc_sum = np.sum(sum_y_zs)
    fill_vol = calc_sum * dx * dy * dz

    summed_porous = (pipe_vol - fill_vol) / pipe_vol
    print("calc from sum:", summed_porous)

    # Calculate the relative error between these methods
    err = abs(summed_porous-porosity)/summed_porous
    print("err:", err)

    rt = time.time() - start_time
    print("--- %s seconds ---" % rt)

    # Now write this data out to a file for later viewing

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




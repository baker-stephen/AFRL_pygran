import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    positions = []
    x_min = 2
    x_max = -2
    y_min = 2
    y_max = -2
    z_min = 4
    z_max = -.2

    with open('legacy/csvs_new/particles5070000.csv', 'r') as r:
        r.readline()
        N_parts = 1000
        for line in r:
            xyz = [float(p) for p in line.split(",")]
            positions.append(xyz)
            if xyz[2] + .5 / 25.4 < 1:
                x_min = min(x_min, xyz[0])
                x_max = max(x_max, xyz[0])
                y_min = min(y_min, xyz[1])
                y_max = max(y_max, xyz[1])
                z_min = min(z_min, xyz[2])
                z_max = max(z_max, xyz[2])
        r.close()

    # print(positions)
    positions = sorted(positions, key=lambda x: x[2])
    # print(positions)
    print("x_min:", x_min)
    print("x_max:", x_max)
    print("y_min:", y_min)
    print("y_max:", y_max)
    print("z_min:", z_min)
    print("z_max:", z_max)

    # sphere of radius .5mm
    ball_r = .5 / 25.4
    ball_r_sq = ball_r**2 #precompute r^2
    full_ball_vol = (4 / 3) * np.pi * ball_r ** 3
    # print("full_ball_vol:", full_ball_vol)
    total_balls_vol = 0

    pipe_in_rad = .26 / 2

    # Only calculate porosity for pellets below 1 in (.9 in fill height), since above that there is irregular shape and settling
    max_z = 1

    x_res = 1000
    x0 = -pipe_in_rad
    xM = pipe_in_rad
    dx = (xM - x0) / (x_res + 1)

    y_res = 1000
    y0 = -pipe_in_rad
    yM = pipe_in_rad
    dy = (yM - y0) / (y_res + 1)

    z_res = 10000
    z0 = .1
    zM = 1.0
    dz = (zM - z0) / (z_res + 1)

    is_filled = []
    for i in range(x_res + 1):
        is_filled.append([])
        for j in range(y_res + 1):
            is_filled[i].append(np.zeros(z_res + 1))

    for p in positions:
        xc = p[0]
        yc = p[1]
        zc = p[2]

        # If the entirety of the current position is greater than height max, we're done, since sorted by z
        if zc - ball_r > max_z:
            break

        # The whole ball is within the region of interest
        elif zc + ball_r < max_z:
            total_balls_vol += full_ball_vol

        # A portion of the ball is within the region of interest, subtract the volume that is above
        else:
            h = zc + ball_r - max_z
            vol_cap = (1 / 3) * np.pi * h ** 2 * (3 * ball_r - h)
            total_balls_vol += full_ball_vol - vol_cap


        #Calculate porosity as a function of x,y,z

        # lowest extent of this ball
        z_low = max(z0, zc - ball_r)
        # nearest logical index to the lowest extent of this ball
        # z_lli = (z_low - z0) / dz
        z_low_li = max(round(-.5 + ((z_low - z0) / dz)), 0)
        # print("diff low: ", z_low_li-z_lli)
        # print("rounded to:",z_low_li)
        # print("rounded from:",z_lli)
        # repeat for highest extent
        z_high = min(zM, zc + ball_r)
        z_high_li = min(round(.5 + ((z_high - z0) / dz)), z_res)
        z_limit = True #True when z distance from center of ball = ball radius
        z_frac = 1.0 #How much of the ball covers this z logical coordinate
        for z in range(z_low_li, z_high_li):
            # given this z coordinate, calculate the x extents:
            z_phys = (z * dz + z0)  # recover physical z
            z_ball = abs((z_phys - z_low) - ball_r)  # logical z distance from center of ball
            z_ball_sq = z_ball ** 2
            x_left_li = 0
            x_right_li = 0
            x_left = xc
            x_right = xc
            if z_ball > ball_r:
                z_frac = 1-((z_ball - ball_r) / dz)
                if z_frac > 1.0:
                    print("z_frac bigger:",z_frac)
                    break
            else:
                z_frac = 1.0
            if z_ball_sq >= ball_r_sq:
                # At the limit of z, lock x here
                # print("percent overlap:", (z_ball-ball_r)/z_ball)
                # print("physical range in radii:", (z_high-z_low)/ball_r)
                # print("z height before sub:",z_phys - z_low)
                # print("ball r:",ball_r)
                # continue
                z_limit = True
                # z_frac = (z_ball-ball_r)/dz
                if z_frac > 1.0:
                    print("z_frac bigger:",z_frac)
                    break
                x_left_li = round((xc - x0) / dx)
                x_right_li = x_left_li+1
            else:
                z_limit = False
                # z_frac = 1.0
                # if z_ball_sq >= ball_r_sq:
                #     print("xr problem")
                #     print("z_ball_sq:", z_ball_sq)
                #     print("ball_r_sq:", ball_r_sq)
                #     print("z_phys:", z_phys)
                #     print("z_ball:", z_ball)
                #     print("z:", z)
                #     print("z_low_li:", z_low_li)
                #     print("z_low:", z_low)
                #     break
                xr = np.sqrt(ball_r_sq - z_ball_sq)
                x_left = max(x0, xc - xr)
                x_left_li = max(round(-.5+((x_left - x0) / dx)), 0)
                x_right = min(xM, xc + xr)
                x_right_li = min(round(.5+((x_right - x0) / dx)), x_res)
            x_frac = 1.0
            for x in range(x_left_li, x_right_li):
                x_phys = (x * dx + x0)
                x_ball = abs((x_phys - x_left) - ball_r)
                x_ball_sq = x_ball**2
                y_back = yc
                y_back_li = 0
                y_front = yc
                y_front_li = 0
                if x_ball > ball_r:
                    x_frac = (x_ball - ball_r) / dx
                    if x_frac > 1.0:
                        print("x_frac bigger:", x_frac)
                        break
                else:
                    x_frac = 1.0
                if z_limit or z_ball_sq + x_ball_sq >= ball_r_sq:
                    # continue
                    # at the x limit, lock y
                    # x_frac = (x_ball - ball_r) / dx
                    # if x_frac >= 1.0:
                    #     print("x_frac bigger")
                    #     break
                    y_back_li = round((yc - y0) / dy)
                    y_back_li = y_back_li+1
                else:
                    # x_frac = 1.0
                    yr = np.sqrt(ball_r_sq - z_ball_sq - x_ball_sq)
                    y_back = max(y0, yc - yr)
                    y_back_li = max(round(-.5+((y_back - y0) / dy)), 0)
                    y_front = min(yM, yc + yr)
                    y_front_li = min(round(.5+((y_front - y0) / dy)), y_res)
                for y in range(y_back_li, y_front_li):
                    y_phys = (y * dy + y0)
                    y_ball = abs((y_phys - y_back) - ball_r)
                    y_frac = 1.0
                    if y_ball > ball_r:
                        y_frac = 1-((y_ball - ball_r) / dy)
                        if y_frac < 0.0:
                            y_frac = 0.0
                            print("y_frac negative:", y_frac)
                            break
                    else:
                        y_frac = 1.0
                    if z_frac*x_frac*y_frac>1.0:
                        print("filling with too big:",z_frac*x_frac*y_frac)
                    elif z_frac*x_frac*y_frac<0:
                        print("filling with negative:", z_frac * x_frac * y_frac)
                    # elif z_frac*x_frac*y_frac!=1.0:
                    #     print("Good! filling with fraction:", z_frac * x_frac * y_frac)
                    is_filled[x][y][z] = z_frac*x_frac*y_frac

    # print("total_balls_vol:", total_balls_vol)

    # print("pipe_in_rad:",pipe_in_rad)
    # print("max x extent:",x_max+ball_r)
    # print("max x extent minus skin:", x_max + ball_r - (.5/25.4)/2)

    pipe_vol = np.pi * pipe_in_rad ** 2 * (max_z - .1)

    # print("OG fill vol:", total_balls_vol)

    porosity = (pipe_vol - total_balls_vol) / pipe_vol

    print("total porosity:", porosity)

    # calculate porosity by averaging over filled to make sure resolutions are adequate
    # avg_zs = []
    # avg_y_zs = []
    # for i in range(x_res):
    #     avg_zs.append([])
    #     for j in range(y_res):
    #         avg_zs[i].append(np.mean(is_filled[i][j]))
    #     avg_y_zs.append(np.mean(avg_zs[i]))
    # calc_avg = np.mean(avg_y_zs)
    # print("calc_avg:",1-calc_avg)

    sum_zs = []
    sum_y_zs = []
    for i in range(x_res+1):
        sum_zs.append([])
        for j in range(y_res+1):
            sum_zs[i].append(np.sum(is_filled[i][j]))
        sum_y_zs.append(np.sum(sum_zs[i]))
    calc_sum = np.sum(sum_y_zs)
    fill_vol = calc_sum * dx * dy * dz

    # domain_vol = 0
    # for x in range(x_res+1):
    #     x_phys = (x * dx + x0)
    #     y_min = -np.sqrt(pipe_in_rad**2-x_phys**2)
    #     y_max = -y_min
    #     y_min_li = max(round(-.5 + ((y_min - y0) / dy)), 0)
    #     y_max_li = min(round(.5 + ((y_max - y0) / dy)), y_res)
    #     for y in range(y_min_li,y_max_li):
    #         domain_vol+=z_res
    #
    # domain_vol *= dx*dy*dz

    # print("fill_vol:", fill_vol)
    print("pipe vol:", pipe_vol)
    # print("domain_vol:", domain_vol)
    summed_porous = (pipe_vol - fill_vol) / pipe_vol
    print("calc from sum:", summed_porous)
    # summed_domain = (domain_vol - fill_vol) / domain_vol
    # print("calc from sum domain:", summed_domain)

    print("err:",abs(summed_porous-porosity)/summed_porous)
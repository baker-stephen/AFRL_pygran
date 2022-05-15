import numpy as np
import matplotlib.pyplot as plt


def xtoli(x: float, x0: float, dx: float) -> float:
    return (x-x0)/dx

if __name__ == "__main__":
    positions = []

    with open('particles_good_final.csv','r') as read:

        #Load final particle positions as polar coords (r, theta, z)
        particle_count = 73
        read.readline().split(',')
        data = read.readline().split(',')
        i = 0
        while data[0] != "-1":
            positions.append([])
            print(data)
            positions[i].append(float(data[0]))
            positions[i].append(float(data[1]))
            positions[i].append(float(data[2]))
            # print(positions[i][2])
            i += 1
            data = read.readline().split(',')
        read.close()

    print("done read")
    #Determine whether a particle is present at each polar coordinate

    #Bounds and detail
    resolution = 100
    zs = np.linspace(.1, 4.1, resolution)
    z0 = zs[0]
    zf = zs[-1]
    dz = (zf-z0)/(resolution+1)
    thetas = np.linspace(0, 2*np.pi, resolution)
    theta0 = thetas[0]
    thetaf = thetas[-1]
    dtheta = (thetaf-theta0)/(resolution+1)
    rs = np.linspace(0, 5/16, resolution)
    r0 = rs[0]
    rf = rs[-1]
    dr = (rf-r0)/(resolution+1)


    #arrays to hold whether a particle is located at a polar position
    filled = []
    for i in range(resolution):
        filled.append([])
        for j in range(resolution):
            filled[i].append(np.zeros(resolution))

    r_ball_sq = (1/8)**2
    r_ball = 1/8
    #Yes I am aware how inefficient this is

    # for j,theta in enumerate(thetas):
    #     costh = np.cos(theta)
    #     sinth = np.sin(theta)
    #     for k,r in enumerate(rs):
    #         x = r * costh
    #         y = r * sinth
    #         for i,z in enumerate(zs):
    #             # convert polar to cartesian
    #             for p in positions:
    #                 d_sq = (x-p[0])**2+(y-p[1])**2+(z-p[2])**2
    #                 if d_sq<=r_ball_sq:
    #                     #there's something here!
    #                     filled[i][j][k] = 1
    #                     break

    for p in positions:
        x_p = p[0]
        y_p = p[1]
        z_p = p[2]
        r_p = np.sqrt(x_p**2+y_p**2)
        r_min = int(xtoli(max(r_p-r_ball, 0.0), r0, dr))
        r_max = int(xtoli(r_p+r_ball, r0, dr)) + 1
        z_min = int(xtoli(max(z_p-r_ball, 0.0), z0, dz))
        z_max = int(xtoli(max(z_p+r_ball, 0.0), z0, dz))+1
        for j, theta in enumerate(thetas):
            costh = np.cos(theta)
            sinth = np.sin(theta)
            k = r_min
            for r in rs[r_min:r_max+1]:
                x = r * costh
                y = r * sinth
                i = z_min
                for z in zs[z_min:z_max+1]:
                    if filled[i][j][k] > .5:
                        continue
                    d_sq = (x - x_p) ** 2 + (y - y_p) ** 2 + (z - z_p) ** 2
                    if d_sq <= r_ball_sq:
                        filled[i][j][k] = 1
                    i += 1
                k += 1

    print("done with determining things")

    avg_zs = []
    avg_theta_zs = []
    for i in range(resolution):
        avg_zs.append([])
        for j in range(resolution):
            avg_zs[i].append(np.mean(filled[i][j]))
        avg_theta_zs.append(np.mean(avg_zs[i]))

    print("print avgs:",avg_theta_zs)

    portions = avg_theta_zs
    porosities = []
    for p in portions:
        p = float(p)
        porosities.append(1.0 - p)

    print(porosities)

    plt.plot(zs, porosities)
    plt.show()
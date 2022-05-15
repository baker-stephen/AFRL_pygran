import numpy as np

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
    resolution = 300
    zs = np.linspace(.1,4.1,resolution)
    thetas = np.linspace(0,2*np.pi,resolution)
    rs = np.linspace(0,5/16,resolution)

    #arrays to hold whether a particle is located at a polar position
    filled = []
    for i in range(resolution):
        filled.append([])
        for j in range(resolution):
            filled[i].append(np.zeros(resolution))

    r_ball_sq = (1/8)**2
    #Yes I am aware how inefficient this is
    for i,z in enumerate(zs):
        for j,theta in enumerate(thetas):
            for k,r in enumerate(rs):
                # convert polar to cartesian
                c = [r * np.cos(theta), r * np.sin(theta), z]
                for p in positions:
                    d_sq = (c[0]-p[0])**2+(c[1]-p[1])**2+(c[2]-p[2])**2
                    if d_sq<=r_ball_sq:
                        #there's something here!
                        filled[i][j][k] = 1

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

    plt.plot(rs, porosities)
    plt.show()
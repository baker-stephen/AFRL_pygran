import numpy as np
import matplotlib.pyplot as plt

if __name__ =="__main__":
    #pre-shake"
    positions = []
    zs = []
    with open('legacy/csvs_new/particles2290000.csv', 'r') as r:
        r.readline()
        for line in r:
            xyz = [float(p) for p in line.split(',')]
            positions.append(xyz)
            zs.append(xyz[2])
        r.close()
    print("zs:",zs)
    # npp = np.array(positions)
    # print("npp:", npp)
    sort_zs = np.sort(zs)
    print("sp:",sort_zs)

    #count number of whole particles below 1":
    Rp = (1/25.4)/2
    print("lowest:",sort_zs[0]-Rp)

    last_index = len(sort_zs)-1
    for i in range(len(sort_zs)-1,0,-1):
        if sort_zs[i]<=1+Rp:
            last_index = i
            break
    whole_parts = sort_zs[:last_index+1]
    print("last_index",last_index)

    vol_sphere = (4/3)*np.pi*Rp**3
    vol_all_spheres = (last_index+1)*vol_sphere
    pipe_vol = .9*np.pi*(.26/2)**2
    poros = (pipe_vol-vol_all_spheres)/pipe_vol
    print("porosity only counting completely contained:")
    print(poros)


    # plt.plot(sort_zs)
    # plt.show()


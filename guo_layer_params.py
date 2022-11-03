import numpy as np
from scipy.optimize import curve_fit
from fluent.fluent import poros_Guo

ID_DP_dict = {'0.26': ['0.7/25.4', '0.8/25.4', '0.9/25.4', '0.85/25.4', '1/25.4', '1.1/25.4', '1.2/25.4', '1.5/25.4', '1.7/25.4'],
              '0.602': ['1/4', '1/8', '1/16', '3/16', '3/32', '7/16'],
              '1.029': ['1/4', '1/8', '3/16', '5/16', '7/16'],
              '1.59': ['1/4', '1/8', '3/16', '5/16', '5/32', '7/16', '7/32', '9/32', '9/64', '15/64']}

if __name__=="__main__":
    zdata = []

    pipe_in_rad = 1.029 / 2

    x_res = 800
    x0 = -pipe_in_rad
    xM = pipe_in_rad
    dx = (xM - x0) / (x_res + 1)

    y_res = 800
    y0 = -pipe_in_rad
    yM = pipe_in_rad
    dy = (yM - y0) / (y_res + 1)

    z_res = 2500
    #TODO: adjust zM,z0
    z0 = .5
    zM = 5.0
    dz = (zM - z0) / (z_res + 1)
    with open('avg_z_'+str(x_res)+'-'+str(y_res)+'-'+str(z_res)+'.npy', 'rb') as znpy:
        zdata = np.load(znpy)
        znpy.close()

    start_phys = 1.0
    start_ind = (start_phys-z0)/dz
    print("start_ind:",start_ind)
    start_ind = round(start_ind)

    minim = 1.0
    maxim = 0.0

    for z in zdata[start_ind:]:
        if z<minim:
            minim = z
        if z>maxim:
            maxim = z

    print("min", minim)
    print("max", maxim)

    amp = (maxim - minim)/2
    print("amp:",amp)
    off = (maxim + minim)/2
    print("off:",off)

    positions = []
    with open('particles.csv','r') as r:
        r.readline()
        for line in r:
            xyz = [float(p) for p in line.split(",")]
            positions.append(xyz)
        r.close()

    positions = sorted(positions, key=lambda x: x[2])
    DP = 1/4
    sphere_rad = DP/2
    ball_r = sphere_rad
    ball_r_sq = ball_r ** 2.0
    full_ball_vol = (4.0 / 3.0) * np.pi * ball_r ** 3.0

    total_balls_vol = 0
    for p in positions:
        xc = p[0]
        yc = p[1]
        zc = p[2]

        # If the entirety of the current position is greater than height max, we're done, since sorted by z
        if zc - ball_r > zM:
            break

        # the entirety of the current ball is below the height minimum, disregard
        elif zc + ball_r < start_phys:
            continue

        # The whole ball is within the region of interest
        elif zc + ball_r < zM and zc - ball_r > start_phys:
            total_balls_vol += full_ball_vol

        # A portion of the ball is within the region of interest, subtract the volume that is above
        elif zc - ball_r > start_phys:
            h = zc + ball_r - zM
            vol_cap = (1.0 / 3.0) * np.pi * h ** 2.0 * (3.0 * ball_r - h)
            total_balls_vol += full_ball_vol - vol_cap

        # A portion of the ball is within the region of interest, subtract the volume that is below
        else:
            h = zc + ball_r - start_phys
            vol_cap = (1.0 / 3.0) * np.pi * h ** 2.0 * (3.0 * ball_r - h)
            total_balls_vol += vol_cap

    pipe_vol = np.pi * pipe_in_rad ** 2 * (zM - start_phys)

    porosity_calc = (pipe_vol - total_balls_vol) / pipe_vol
    print("calculated porosity:",porosity_calc)

    fguesses = []
    errs = []
    phis = []
    lhs = []
    k1s = []
    k2s = []
    mzcs = []
    for guess in range(45,55):
        # guess = 50
        def to_fit(x: float, f: float, phi: float):
            return amp*np.sin(guess*2*np.pi*f*x/(zM-start_phys)+phi)+off


        xdata = np.linspace(start_phys,zM,len(zdata[start_ind:]))
        print("xdata shape:", xdata.shape)
        print("xdata:", xdata)
        ydata = np.array(zdata[start_ind:])
        print("ydata shape:", ydata.shape)
        print(ydata)

        coefs, something = curve_fit(to_fit, xdata, ydata,maxfev=10000)
        print("coefs:", coefs)
        print("f:",coefs[0])
        print("f guess:", coefs[0]*guess)
        fguesses.append(coefs[0]*guess)
        print("phi:",coefs[1])
        phis.append(coefs[1])



        layer_height = 2*(zM-start_phys)/(guess*coefs[0])
        print("layer h:",layer_height)
        lhs.append(layer_height)
        print("sph rad:",sphere_rad)
        k2 = DP/layer_height
        print("k2:",k2)
        k2s.append(k2)



        z_counts = np.zeros(len(zdata[start_ind:]))
        for i,z_ind in enumerate(range(len(zdata[start_ind:]))):
            z = z0 + (i + start_ind) * dz
            for pos in positions:
                if z-(layer_height/2) < pos[2] < z+(layer_height/2):
                    z_counts[z_ind]+=1
        print("zcounts:",z_counts)
        mzc = np.mean(z_counts)
        mzcs.append(mzc)
        print("mean zcounts:",mzc)

        k1 = mzc*sphere_rad**2/(pipe_in_rad**2)
        k1s.append(k1)
        print("k1:",k1)

        poros = 1 - (2/3)*k1*k2
        print("poros:",poros)
        errs.append(abs(poros-porosity_calc)/porosity_calc)
        print("err:",abs(poros-porosity_calc)/porosity_calc)
    N = pipe_in_rad / sphere_rad
    print("N:", N)
    print("Guo poros:",poros_Guo(N))

    print("errs")
    print(errs)
    min_err = 1.0
    err_ind = -1
    for i,err in enumerate(errs):
        if err<min_err:
            min_err = err
            err_ind = i

    print("min_err:",min_err)
    print("best fguess:",fguesses[err_ind])
    print("best phi:",phis[err_ind])
    print("best mzc:", mzcs[err_ind])
    print("best lh:", lhs[err_ind])


    with open('guo_layers.csv','w') as csv:
        csv.write('z,data,fit\n')
        for i,z in enumerate(zdata[start_ind:]):
            csv.write(str(z0+(i+start_ind)*dz)+',')
            csv.write(str(z)+',')
            csv.write(str(to_fit(z0+(i+start_ind)*dz,fguesses[err_ind],phis[err_ind]))+'\n')

        csv.close()
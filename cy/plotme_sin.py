import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import statistics as stat
import time
from statsmodels.regression.linear_model import OLS

import sys
sys.path.append('../')
from param_defn import PD

def fluent_print(params: list):
    ret = str(params[0])
    for i in range(1, len(params)):
        if i % 2 != 0:
            ret += "+"+str(params[i])+"*cos(2[in^-1]*PI*"+str(i // 2 + 1) + "*p_r)"
        else:
            ret += "+" + str(params[i]) + "*sin(2[in^-1]*PI*" + str(i // 2) + "*p_r)"
    return ret

def sin_cos_list(x,params):
    ret = params[0]
    for i in range(1,len(params)):
        if i % 2 != 0:
            ret += params[i]*np.cos(2 * np.pi * (i//2 + 1) * x)
        else:
            ret += params[i]*np.sin(2 * np.pi * (i//2) * x)
    return ret

def go(input_params: PD):

    wrt_dir = input_params.out_dir
    z_0 = input_params.z0()
    z_M = input_params.zM()
    res = input_params.res()
    ID = input_params.ID


    print("start")
    start_time = time.time()
    out = open(wrt_dir+'outputs.txt',  'a')
    # radial = []
    pipe_in_rad = ID / 2
    x_res = res[0]
    x0 = -pipe_in_rad
    xM = pipe_in_rad
    dx = (xM - x0) / (x_res + 1)
    y_res = res[1]
    y0 = -pipe_in_rad
    yM = pipe_in_rad
    dy = (yM - y0) / (y_res + 1)
    z_res = res[2]
    z0 = z_0
    zM = z_M
    dz = (zM - z0) / (z_res + 1)
    r_res = res[3]
    rs = np.linspace(0, pipe_in_rad, r_res+1)
    print("opening avg r arr")
    with open(wrt_dir+'avg_r_' + str(x_res) + '-' + str(y_res) + '-' + str(z_res) + '-' + str(r_res) + '.npy', 'rb') as f:
        radial = np.load(f)
        f.close()
    print("done")
    print("std dev r:",stat.stdev(radial[1:]))
    out.write("std dev r, after first index: " + str(stat.stdev(radial[1:])) + '\n')

    start_ind = 1
    for i,r in enumerate(radial):
        if 0>=r or r>=1:
            print("bad poros:",r)
            start_ind=i+1
    print("start_ind:",start_ind)
    out.write("start_ind: " + str(start_ind) + '\n')
    print("start poros:",radial[start_ind])
    out.write("start poros: " + str(radial[start_ind]) + '\n')
    print("start radius:",rs[start_ind])
    out.write("start radius: " + str(rs[start_ind]) + '\n')
    for i in range(start_ind):
        radial[i]=radial[start_ind]
    fig, ax = plt.subplots()
    ax.set_xlabel("radius (in)")
    # dp = Dp
    ax.plot(rs, radial,label='data')

    ny = len(radial)
    N = 40
    nparams = 2*N + 1

    xas = []
    for i in range(1,N+1):
        xas.append(np.cos(2 * np.pi * i * rs))

    xbs = []
    for i in range(1,N+1):
        xbs.append(np.sin(2 * np.pi * i * rs))

    matr = np.ones((ny, nparams), dtype=float)
    for i in range(1,nparams):
        if i%2!=0:
            matr[:,i] = xas[i//2]
        else:
            matr[:, i] = xbs[i // 2 - 1]

    model = OLS(radial, matr)
    results = model.fit()
    print(results.params)
    print(results.rsquared)
    print(results.resid)
    out.write("Fourier fit degree " + str(N) + '-\n')
    out.write(fluent_print(results.params) + '\n')
    out.write("Residuals abs avg: " + str(np.mean(abs(results.resid))) + '\n')
    out.write("R squared: " + str(results.rsquared) + '\n')

    ax.plot(rs, sin_cos_list(rs, results.params), 'r-',
            label='Fourier fit deg '+str(N))
    ax.legend()
    plt.savefig(wrt_dir+"r-poros-fit_" + str(x_res) + '-' + str(y_res) + '-' + str(z_res) + '-' + str(r_res) + '.png')
    # plt.show()



    def f1(x):
        return sin_cos_list(x,results.params)
    fit_sin_cos_integ = quad(f1, 0.0, pipe_in_rad)
    print("fit_one_integ:", fit_sin_cos_integ)
    print("avg poros fit one:", fit_sin_cos_integ[0] / pipe_in_rad)
    out.write("avg poros fit integ: " + str(fit_sin_cos_integ[0] / pipe_in_rad) + '\n')
    n = 10000
    r_fits = np.linspace(0,pipe_in_rad,n)
    p_fits = []
    for i,r in enumerate(r_fits):
        p_fits.append(sin_cos_list(r, results.params))
    print("mean from generated pts:",np.mean(p_fits))
    out.write("mean from generated pts: " + str(np.mean(p_fits)) + '\n')
    print("std dev from generated pts:",np.std(p_fits))
    out.write("std dev from generated pts: " + str(np.std(p_fits)) + '\n')
    rt = time.time()-start_time
    print("--- %s seconds ---" % rt)
    out.write('plotme runtime: ' + str(rt) + '\n')

    def f1(x):
        return sin_cos_list(x, results.params) * x

    # vol avg poros = 2*pi*L*integral 0->R of poros(r)*r*dr / pi*R**2*L
    fit_sin_cos_integ = quad(f1, 0.0, float(ID) / 2)
    print("integral:", fit_sin_cos_integ)
    vol_avg_poros = 2 * fit_sin_cos_integ[0] / ((float(ID) / 2) ** 2)
    print("volume averaged porosity:", vol_avg_poros)
    out.write("\nvolume averaged porosity: " + str(vol_avg_poros))
    out.close()
    print("Done!")
from dimension import *
import matplotlib.pyplot as plt

from scipy.optimize import fsolve


def tau_ann(l: Length, d: Length) -> float:
    """
    Tortuosity in the annulus. Calculated by averaging the path length when traveling directly upwards until meeting
    the sphere, then flowing along sphere for varying initial x positions.
    :param l: layer overlap (y difference between sphere centers in consecutive layers)
    :param d: sphere diameter
    :return: tau_a
    """
    # return float(1 + d * (np.pi ** 2 - 8) / (8 * np.pi * l))
    return float(1 + 0.5*d*(1-np.pi/4)/l)


def tau_core(l: Length, d: Length) -> float:
    """
    Tortuosity in the core. Calculated by calculating the path length when traveling directly upwards until meeting
    the sphere, or starting flow along sphere and then flowing directly upwards.
    :param l: layer overlap (y difference between sphere centers in consecutive layers)
    :param d: sphere diameter
    :return: tau_c
    """
    return float(1 + 0.5 * d * (np.pi / 2 - 1) / l)


def A_intercept(D: Length, d: Length) -> Derived:
    """
    The area of the circle which intersects a larger circle
    :param d: sphere diameter (m)
    :param D: Inner diameter of the tube
    :return: A_int (m^2)
    """
    # print("D:",D)
    # print("Dp:",d)
    r = d / 2
    R = D / 2 - r
    d1 = (2 * R ** 2 - r ** 2) / (2 * R)
    d2 = R - d1
    A_int = R ** 2 * np.arccos(d1 / R) - d1 * np.sqrt(R ** 2 - d1 ** 2) + r ** 2 * np.arccos(d2 / r) \
            - d2 * np.sqrt(r ** 2 - d2 ** 2)
    return A_int


def Dh_core(D: Length, d: Length, n: float, na: float) -> Length:
    """
    Mean capillary tube diameter in the core.
    :param D: Inner diameter of the tube
    :param d: Sphere diameter
    :param n: Number of spheres per layer
    :param na: Number of spheres per layer forming the annulus
    :return: D_h_c
    """
    N = float(D / d)
    A_int = A_intercept(D, d)
    return d * np.sqrt((N ** 2 + 2 * N - n + na + 1 - 4 * na * A_int / (np.pi * d ** 2)) / (n-1))

def Acore(D: Length, d: Length, n: float, na: float) -> Derived:
    """
    Mean capillary tube area in the core.
    :param D: Inner diameter of the tube
    :param d: Sphere diameter
    :param n: Number of spheres per layer
    :param na: Number of spheres per layer forming the annulus
    :return: D_h_c
    """
    N = float(D / d)
    A_int = A_intercept(D, d)
    return 0.25 * np.pi * d ** 2 * (N ** 2 + 2 * N - n + na + 1 - 4 * na * A_int / (np.pi * d ** 2)) / (n - 1)

def Atotal_core(D: Length, d: Length, n: float, na: float) -> Derived:
    """
    Mean capillary tube area in the core.
    :param D: Inner diameter of the tube
    :param d: Sphere diameter
    :param n: Number of spheres per layer
    :param na: Number of spheres per layer forming the annulus
    :return: D_h_c
    """
    N = float(D / d)
    A_int = A_intercept(D, d)
    return 0.25 * np.pi * d ** 2 * (N ** 2 + 2 * N - n + na + 1) - na * A_int

def Apath_core(d: Length)->Derived:
    R = d/2
    return R**2*(4-np.pi)

def Pcore(Acore: Derived, d: Length)->Length:
    """
    Capillary tube perimeter in the core
    :param Acore: Mean capillary tube area in the core (m^2)
    :param d: Sphere diameter
    :return: Pc
    """
    R = float(d)/2
    def f(P):
        x0 = R*np.cos(P[0]/(4*R))
        Acirc_int = np.pi*R**2/4 - 0.5*x0*np.sqrt(R**2-x0**2) - 0.5*R**2*np.arctan(x0/np.sqrt(R**2-x0**2))
        s = 2*(x0-2*R/np.sqrt(5))
        A_calc = s**2 - 8*Acirc_int
        return float(Acore)-A_calc

    root = fsolve(f,float(np.sqrt(Acore)))
    print("root:",root)
    print(f(root))
    print(np.isclose(f(root), 0.0))


def Dh_ann(D: Length, d: Length, na: float) -> Length:
    """
    Mean capillary tube diameter in the annulus.
    :param D: Inner diameter of the tube
    :param d: Sphere diameter
    :param na: Number of spheres per layer forming the annulus
    :return: D_h_a
    """
    A_int = A_intercept(D, d)
    # return np.sqrt((2*d)/na * (D - d/2) - d**2 + 4*A_int/np.pi)
    return d * np.sqrt(1 / na * (2 * D / d - 1) - 1 + 4 * A_int / (np.pi * d ** 2))


def alpha(cf1: float, cf2: float, rho: Derived, tau: float, mu: Derived, u: Derived, D_H: Length,
          phi: float) -> Derived:
    """
    Viscous/skin resistance coefficient, appears next to viscosity times superficial velocity.
    :param cf1: Scaling coefficient for skin friction coefficient versus Reynolds number
    :param cf2: Exponential coefficient for skin friction coefficient versus Reynolds number
    :param rho: mass density (kg/m^3)
    :param tau: tortuosity
    :param mu: dynamic viscosity (kg/m*s)
    :param u: superficial velocity (m/s)
    :param D_H: Mean flow path diameter (m)
    :param phi: porosity
    :return: alpha (1/m^2)
    """
    return 2 * cf1 * rho ** (1 - cf2) * tau ** (2 - cf2) * mu ** (cf2 - 1) * u ** (1 - cf2) / (
            D_H ** (1 + cf2) * phi ** (2 - cf2))

def alpha_new(cf1: float, cf2: float, rho: Derived, tau: float, mu: Derived, u: Derived, d: Length,
          phi: float, a: float) -> Derived:
    """
    Viscous/skin resistance coefficient, appears next to viscosity times superficial velocity.
    :param cf1: Scaling coefficient for skin friction coefficient versus Reynolds number
    :param cf2: Exponential coefficient for skin friction coefficient versus Reynolds number
    :param rho: mass density (kg/m^3)
    :param tau: tortuosity
    :param mu: dynamic viscosity (kg/m*s)
    :param u: superficial velocity (m/s)
    :param d: particle diameter (m)
    :param phi: porosity
    :param a: area expansion factor
    :return: alpha (1/m^2)
    """
    hydr_D = 0.5*d*(4*a**2-np.pi)/(2*(a-1) + np.pi/2)
    return cf1 * 2 * np.pi * rho ** (1 - cf2) * tau ** (2 - cf2) * mu ** (cf2 - 1) * u ** (1 - cf2) / (
        2 * (a**2 - np.pi/4) * phi ** (2 - cf2) * d * hydr_D**cf2)

def alpha_ann(cf1: float, cf2: float, rho: Derived, tau: float, mu: Derived, u: Derived, d: Length,
          phi: float, D: Length) -> Derived:
    Aint = A_intercept(D,d)
    N = float(D/d)
    theta = theta_ann(N)
    n_a = na(N)
    Theta_ex = (2/n_a)*(np.pi - n_a*theta)
    theta_adj = theta + 0.5*Theta_ex
    hydr_D = 4*(d**2*theta_adj*(2*N-1)/4 - (np.pi*0.25*d**2 - Aint))/(theta_adj*D + 0.5*np.pi*d + Theta_ex*(D/2 - d/2))
    return cf1 * (theta_adj*D + 0.5*np.pi*d) * rho ** (1 - cf2) * tau ** (2 - cf2) * mu ** (cf2 - 1) * u ** (1 - cf2) / (
        2 * (d**2*theta_adj*(2*N-1)/4 - (np.pi*0.25*d**2 - Aint)) * phi ** (2 - cf2) * hydr_D**cf2)

def expand_factor_core(d: Length, n: float, D: Length, n_a: int)->float:
    """
    Area expansion factor
    :param d:
    :param n:
    :param D:
    :param n_a:
    :return: a
    """
    At = Atotal_core(D,d,n,n_a)
    N_guess = 1
    if n>3:
        N_guess = n - 1
    Aguess = At/N_guess
    a = np.sqrt((Aguess+np.pi*d**2/4)/d**2)
    return float(a)

def beta(cD: float, tau: float, D_H: Length, phi: float) -> Derived:
    """
    Inertial/orifice resistance coefficient, appears next to density times superficial velocity squared.
    :param cD: Orifice flow resistance coefficient
    :param tau: tortuosity
    :param D_H: Mean flow path diameter (m)
    :param phi: porosity
    :return: beta (1/m)
    """
    return cD * tau ** 2 / (2 * phi ** 2 * D_H)

def beta_new(cD: float, tau: float, l: Length, phi: float) -> Derived:
    """
    Inertial/orifice resistance coefficient, appears next to density times superficial velocity squared.
    :param cD: Orifice flow resistance coefficient
    :param tau: tortuosity
    :param l: Layer height (m)
    :param phi: porosity
    :return: beta (1/m)
    """
    return cD * tau**2 / (2 * phi ** 2 * l)

def Ergun_alpha(A_E: float, epsilon: float, d: Length) -> Derived:
    """
    Viscous/skin resistance coefficient, appears next to viscosity times superficial velocity.
    :param A_E: dimensionless constant
    :param epsilon: porosity
    :param d: sphere diameter (m)
    :return: Ergun alpha (1/m^2)
    """
    return A_E * (1 - epsilon) ** 2 / (epsilon ** 3 * d ** 2)


def Ergun_beta(B_E: float, epsilon: float, d: Length) -> Derived:
    """
    Inertial resistance coefficient, appears next to density times superficial velocity squared.
    :param B_E: dimensionless constant
    :param epsilon: porosity
    :param d: sphere diameter (m)
    :return: Ergun beta (1/m)
    """
    return B_E * (1 - epsilon) / (epsilon ** 3 * d)

def theta_ann(N: float)->float:
    return np.arccos(1 / np.sqrt((N ** 2 - 2 * N + 1) / (N ** 2 - 2 * N)))

def na(N: float) -> int:
    if N <= 2:
        return 1
    else:
        return int(np.pi / theta_ann(N))


def na_float(N: float) -> float:
    if N <= 2:
        return 1
    else:
        return np.pi / theta_ann(N)


def Guo_params(N: float, d: Length) -> (float, float):
    """
    Guo determined layer height and particles per layer
    :param N: D/d ratio
    :return: (n, l)
    """
    if N >= 3.0 or 1.866 < N < 2.0:
        print("out of range for Guo")
        return -1

    k1 = 1
    k2 = 1

    if N <= 1.866:
        k1 = N ** (-2)
        k2 = 1 / (np.sqrt(1 - (N - 1) ** 2))
    elif 2 <= N <= 3:
        n = 2
        alpha = np.pi / 4
        if 2.1547 <= N < 2.4142:
            n = 3
            alpha = np.pi / 6
        elif 2.4142 <= N < 2.7013:
            n = 4
            alpha = 22.5 * (np.pi / 180)
        elif 2.7013 <= N:
            n = 5
            alpha = 18 * (np.pi / 180)
        k1 = n / (N ** 2)
        # Guo says there should be a Dp in there, but units don't work??
        k2 = 1 / np.sqrt(1 - ((N - 1) * np.sin(alpha)) ** 2)

    return k1 * N ** 2, d / k2


def poros_Guo(N):
    # Analytically derived. Should be exact for perfectly settled/ordered spheres. Very limited range of applicability.
    if N >= 3.0 or 1.866 < N < 2.0:
        # print("out of range for Guo")
        return -1

    k1 = 1
    k2 = 1

    if N <= 1.866:
        k1 = N ** (-2)
        k2 = 1 / (np.sqrt(1 - (N - 1) ** 2))
    elif 2 <= N <= 3:
        n = 2
        alpha = np.pi / 4
        if 2.1547 <= N < 2.4142:
            n = 3
            alpha = np.pi / 6
        elif 2.4142 <= N < 2.7013:
            n = 4
            alpha = 22.5 * (np.pi / 180)
        elif 2.7013 <= N:
            n = 5
            alpha = 18 * (np.pi / 180)
        k1 = n / (N ** 2)
        # Guo says there should be a Dp in there, but units don't work??
        k2 = 1 / np.sqrt(1 - ((N - 1) * np.sin(alpha)) ** 2)

    return 1 - (2 / 3) * k1 * k2

def A_E_B_E(porosity,N,a1,a2,b1,b2,m):
    A_E = (a1 + (a2*porosity/(1-porosity))*(N/(N-1))**m)
    B_E = (b1*np.power((1-porosity)/porosity,1/3)+b2*(N/(N-1))**m)
    return A_E,B_E

def Aw_derive(alpha:Derived,epsilon:float,d:Length,N:float)->float:
    M = M_factor(N,epsilon)
    return float(alpha*epsilon**3*d**2/(M**2*(1-epsilon)**2))

def Bw_derive(beta:Derived,epsilon:float,d:Length,N:float)->float:
    M = M_factor(N,epsilon)
    return float(beta*epsilon**3*d/(M*(1-epsilon)))

def M_factor(N:float,epsilon:float)->float:
    return 1 + 2 / (3 * N * (1 - epsilon))

def reynold(rho: Derived, u: Derived, tau: float, porosity: float, d: Length, mu: Derived, D:Length, is_core: bool, extra: float)->float:
    hydr_D = 0
    if is_core:
        hydr_D = 0.5 * d * (4 * extra ** 2 - np.pi) / (2 * (extra - 1) + np.pi / 2) #core
        # print("core hydr:", hydr_D)
    else:
        Aint = A_intercept(D,d)
        N = float(D/d)
        n_a = na(N)
        Theta_ex = (2 / n_a) * (np.pi - n_a * extra)
        # print("Theta_ex:",Theta_ex*180/np.pi)
        theta_adj = extra + 0.5 * Theta_ex
        hydr_D = 4*(d ** 2 * theta_adj * (2 * N - 1) / 4 - (np.pi * 0.25 * d ** 2 - Aint)) / (
                    theta_adj * D + 0.5 * np.pi * d + Theta_ex * (D / 2 - d / 2))
        # print("annulus hydr:",hydr_D)
    vx = tau*u/porosity
    return rho*vx*hydr_D/mu

# TODO: Differing porosity calculation at core and annulus
# TODO: define annulus hydr

if __name__ == "__main__":
    Ds = []
    ds = []
    Ns = []
    poroses = []
    ns = []
    ls = []
    all_params = []
    with open("../guo_layers.csv", 'r') as csv:
        csv.readline()
        for line in csv:
            params = [float(item.strip()) for item in line.split(',')]
            if params[2]<=2: #N<2, error in Aint calc
                continue
            all_params.append(params)

        csv.close()

    all_params = sorted(all_params, key=lambda x: x[2])
    for params in all_params:
        Ds.append(Length(params[0], Length.inch))
        ds.append(Length(params[1], Length.inch))
        Ns.append(params[2])
        poroses.append(params[3])
        ns.append(params[4])
        ls.append(Length(params[6], Length.inch))


    tau_as = [tau_ann(l, d) for l, d in zip(ls, ds)]
    tau_cs = [tau_core(l, d) for l, d in zip(ls, ds)]

    n_as = [na(N) for N in Ns]

    rho = Derived(998.2, units={Unit(Mass, Mass.kg): 1, Unit(Length, Length.m): -3})
    mu = Derived(0.001003,
                 units={Unit(Mass, Mass.kg): 1, Unit(Length, Length.m): -1, Unit(Time, Time.s): -1})

    vol = Volume(1, Volume.gal)
    duration = Time(1, Time.minute)
    mdot = rho * vol / duration
    areas = [np.pi * D ** 2 / 4 for D in Ds]
    us = [mdot / (rho * area) for area in areas]

    #From wikipedia...
    cf1 = 32
    cf2 = 1
    cD = 0.6+0.85 #avg discharge coef * 2

    # A_E = 150
    # B_E = 1.75
    a1 = 185
    a2 = 17
    b1 = 1.3
    b2 = .03
    m = 2
    aebes = [A_E_B_E(poros,N,a1,a2,b1,b2,m) for poros,N in zip(poroses, Ns)]
    A_Es = [aebe[0] for aebe in aebes]
    B_Es = [aebe[1] for aebe in aebes]

    expand_cores = [expand_factor_core(d,n,D,n_a) for d,n,D,n_a in zip(ds,ns,Ds,n_as)]
    # print("expansion factors:",expand_cores)

    Res_core = [reynold(rho, u, tau, por, d, mu, D, True, a) for u,tau,por,d,D,a in zip(us, tau_cs, poroses, ds, Ds, expand_cores)]
    Res_ann = [reynold(rho, u, tau, por, d, mu, D, False, theta_ann(float(D/d))) for u, tau, por, d, D in zip(us, tau_as, poroses, ds, Ds)]
    # [print("Reynold core:",rec,"annulus:",rea) for rec,rea in zip(Res_core,Res_ann)]

    alpha_cores = [alpha_new(cf1, cf2, rho, tau, mu, u, d, poros, a) for tau, u, d, poros, a in zip(tau_cs, us, ds, poroses, expand_cores)]
    # alpha_cores = [alpha(cf1, cf2, rho, tau, mu, u, Dh, poros) for tau, Dh, poros, a in zip(tau_cs, Dh_cs, poroses, expand_cores)]
    Aw_cores = [Aw_derive(alph,epsilon,d,N) for alph,epsilon,d,N in zip(alpha_cores,poroses,ds,Ns)]
    alpha_anns = [alpha_ann(cf1, cf2, rho, tau, mu, u, d, poros, D) for tau, u, d, poros,D in zip(tau_as, us, ds, poroses,Ds)]
    # alpha_anns = [alpha(cf1, cf2, rho, tau, mu, u, Dh, poros) for tau, Dh, poros, a in zip(tau_as, Dh_as, poroses, expand_cores)]
    Aw_anns = [Aw_derive(alph,epsilon,d,N) for alph,epsilon,d,N in zip(alpha_anns,poroses,ds,Ns)]

    Ergun_alphas = [Ergun_alpha(A_E, poros, d) for A_E, poros, d in zip(A_Es, poroses, ds)]

    beta_cores = [beta_new(cD, tau, l, poros) for tau, l, poros in zip(tau_cs, ls, poroses)]
    Bw_cores = [Bw_derive(bet,epsilon,d,N) for bet,epsilon,d,N in zip(beta_cores,poroses,ds,Ns)]
    beta_anns = [beta_new(cD, tau, l, poros) for tau, l, poros in zip(tau_as, ls, poroses)]
    Bw_anns = [Bw_derive(bet,epsilon,d,N) for bet,epsilon,d,N in zip(beta_anns,poroses,ds,Ns)]

    Ergun_betas = [Ergun_beta(B_E, poros, d) for B_E, poros, d in zip(B_Es, poroses, ds)]

    ann_fs = [2/N - 1/N**2 for N in Ns]
    weight_aws = [af*awa + (1-af)*awc for af,awa,awc in zip(ann_fs, Aw_anns, Aw_cores)]
    weight_bws = [af*bwa + (1-af)*bwc for af,bwa,bwc in zip(ann_fs, Bw_anns, Bw_cores)]


    # fig1, ax1 = plt.subplots()
    # Cheng_Aws = [A_E/M_factor(N,poros)**2 for A_E,N,poros in zip(A_Es,Ns,poroses)]
    # ax1.scatter(Ns,weight_aws,label="Flow paths (weight avg)")
    # ax1.scatter(Ns, Aw_cores, label="Core")
    # ax1.scatter(Ns, Aw_anns, label="Annulus")
    # ax1.plot(Ns,Cheng_Aws,label="Cheng",c='r')
    # ax1.set_xlabel("D/d ratio")
    # ax1.set_ylabel("Aw")
    # ax1.legend()
    # plt.show()

    fig2, ax2 = plt.subplots()
    Cheng_Bws = [B_E / M_factor(N, poros) for B_E, N, poros in zip(B_Es, Ns, poroses)]
    ax2.scatter(Ns,weight_bws,label="Flow paths (weight avg)")
    ax2.scatter(Ns, Bw_cores, label="Core")
    ax2.scatter(Ns, Bw_anns, label="Annulus")
    ax2.plot(Ns, Cheng_Bws,label="Cheng",c='r')
    ax2.set_xlabel("D/d ratio")
    ax2.set_ylabel("Bw")
    ax2.legend()
    plt.show()

    # fig3, ax3 = plt.subplots()
    # ax3.scatter(Ns, ns, label="number per layer")
    # ax3.set_xlabel("D/d ratio")
    # ax3.set_ylabel("quantity")
    # ax3.legend()
    # plt.show()

    # fig4, ax4 = plt.subplots()
    # ax4.scatter(Ns, ls, label="layer heights")
    # ax4.set_xlabel("D/d ratio")
    # ax4.set_ylabel("quantity")
    # ax4.legend()
    # plt.show()

    # fig5, ax5 = plt.subplots()
    # ax5.plot(Ns, tau_as, label="annular tortuosity ")
    # ax5.plot(Ns, tau_cs,label="core tortuosity")
    # ax5.set_xlabel("D/d ratio")
    # ax5.set_ylabel("tortuosity")
    # ax5.legend()
    # plt.show()

    # fig6, ax6 = plt.subplots()
    # ax6.scatter(Ns, [n*(d**2)/(D**2) for n,d,D in zip(ns,ds,Ds)], label="k1")
    # ax6.set_xlabel("D/d ratio")
    # ax6.set_ylabel("quantity")
    # ax6.legend()
    # plt.show()
    #
    # fig7, ax7 = plt.subplots()
    # ax7.scatter(Ns, [d/l for l,d in zip(ls,ds)], label="k2")
    # ax7.set_xlabel("D/d ratio")
    # ax7.set_ylabel("quantity")
    # ax7.legend()
    # plt.show()

    # k1s_pg = []
    # k2s_pg = []
    # Ns_pg = []
    # ds_pg = []
    # Ds_pg = []
    #
    # for i,N in enumerate(Ns):
    #     if 2<N<3:
    #         Ns_pg.append(N)
    #         ds_pg.append(ds[i])
    #         Ds_pg.append(Ds[i])
    #         k1s_pg.append(ns[i]*(ds[i]**2)/(Ds[i]**2))
    #         k2s_pg.append(ds[i]/ls[i])
    #
    # ns_ls_guo = [Guo_params(N, d) for N, d in zip(Ns_pg, ds_pg)]
    # ns_guo = [nl[0] for nl in ns_ls_guo]
    # ls_guo = [nl[1] for nl in ns_ls_guo]
    #
    # fig8, ax8 = plt.subplots()
    # ax8.scatter(Ns_pg, k1s_pg, label="k1s pg")
    # ax8.scatter(Ns_pg, [n*(d**2)/(D**2) for n,d,D in zip(ns_guo,ds_pg,Ds_pg)], label="k1s guo")
    # ax8.set_xlabel("D/d ratio")
    # ax8.set_ylabel("quantity")
    # ax8.legend()
    # plt.show()
    #
    # fig9, ax9 = plt.subplots()
    # ax9.scatter(Ns_pg, k2s_pg, label="k2s pg")
    # ax9.scatter(Ns_pg, [d/l for l,d in zip(ls_guo,ds_pg)], label="k2s guo")
    # ax9.set_xlabel("D/d ratio")
    # ax9.set_ylabel("quantity")
    # ax9.legend()
    # plt.show()

    Ns_Ds = {0.26:[],0.602:[],1.029:[],1.59:[]}
    beta_Ds = {0.26:[],0.602:[],1.029:[],1.59:[]}
    bw_Ds = {0.26:[],0.602:[],1.029:[],1.59:[]}
    tau_Ds = {0.26:[],0.602:[],1.029:[],1.59:[]}
    l_Ds = {0.26: [], 0.602: [], 1.029: [], 1.59: []}
    # print(Ds)
    for i,N in enumerate(Ns):
        Ns_Ds[float(Ds[i].get(Length.inch)).__round__(3)].append(N)
        beta_Ds[float(Ds[i].get(Length.inch)).__round__(3)].append(beta_cores[i])
        bw_Ds[float(Ds[i].get(Length.inch)).__round__(3)].append(Bw_cores[i])
        tau_Ds[float(Ds[i].get(Length.inch)).__round__(3)].append(tau_cs[i])
        l_Ds[float(Ds[i].get(Length.inch)).__round__(3)].append(ls[i])

    fig10, ax10 = plt.subplots()
    for D in Ns_Ds.keys():
        # ax10.scatter(Ns_Ds[D], beta_Ds[D], label="beta core, D="+str(D))
        ax10.scatter(Ns_Ds[D], l_Ds[D], label="l, D=" + str(D))
    # ax10.scatter(Ns, beta_anns, label="beta ann")
    ax10.set_xlabel("D/d ratio")
    ax10.set_ylabel("beta")
    ax10.legend()
    plt.show()


if __name__ == "__ain__":
    D = Length(1, Length.inch)
    N_calcs = 100
    Ns = np.linspace(2.05, 2.95, N_calcs)
    ds = [D / N for N in Ns]
    # print("ds:",ds)
    ns_ls = [Guo_params(N, d) for N, d in zip(Ns, ds)]
    print(ns_ls)
    ns = [nl[0] for nl in ns_ls]
    print(ns)
    ls = [nl[1] for nl in ns_ls]
    # [print(2 * l / d) for l, d in zip(ls, ds)]
    tau_as = [tau_ann(l, d) for l, d in zip(ls, ds)]
    tau_cs = [tau_core(l, d) for l, d in zip(ls, ds)]
    n_as = [na(N) for N in Ns]
    Dh_as = [Dh_ann(D, d, n_a) for d, n_a in zip(ds, n_as)]
    print("dh ann:",Dh_as)
    Dh_cs = [Dh_core(D, d, n, n_a) for d, n, n_a in zip(ds, ns, n_as)]
    print("dh core:", Dh_cs)
    poroses = [poros_Guo(N) for N in Ns]

    # Acores = [Acore(D,d,n,na) for d,n,na in zip(ds,ns,n_as)]
    # print("Acores:",Acores)
    # Pcores = [Pcore(Ac,d) for Ac,d in zip(Acores,ds)]
    # print("Pcores:",Pcores)

    rho = Derived(998.2, units={Unit(Mass, Mass.kg): 1, Unit(Length, Length.m): -3})
    mu = Derived(0.001003,
                 units={Unit(Mass, Mass.kg): 1, Unit(Length, Length.m): -1, Unit(Time, Time.s): -1})

    vol = Volume(1, Volume.gal)
    duration = Time(1, Time.minute)
    mdot = rho * vol / duration
    area = np.pi * D ** 2 / 4
    u = mdot / (rho * area)
    #From wikipedia...
    cf1 = 0.664
    cf2 = 0.5
    cD = 0.6+0.85 #avg discharge coef * 2

    # A_E = 150
    # B_E = 1.75
    a1 = 185
    a2 = 17
    b1 = 1.3
    b2 = .03
    m = 2
    aebes = [A_E_B_E(poros,N,a1,a2,b1,b2,m) for poros,N in zip(poroses,Ns)]
    A_Es = [aebe[0] for aebe in aebes]
    B_Es = [aebe[1] for aebe in aebes]

    expand_cores = [expand_factor_core(d,n,D,n_a) for d,n,n_a in zip(ds,ns,n_as)]
    print("expansion factors:",expand_cores)

    Res_core = [reynold(rho, u, tau, por, d, mu, D, True, a) for tau,por,d,a in zip(tau_cs, poroses, ds, expand_cores)]
    Res_ann = [reynold(rho, u, tau, por, d, mu, D, False, theta_ann(D/d)) for tau, por, d in zip(tau_as, poroses, ds)]
    [print("Reynold core:",rec,"annulus:",rea) for rec,rea in zip(Res_core,Res_ann)]

    alpha_cores = [alpha_new(cf1, cf2, rho, tau, mu, u, d, poros, a) for tau, d, poros, a in zip(tau_cs, ds, poroses, expand_cores)]
    # alpha_cores = [alpha(cf1, cf2, rho, tau, mu, u, Dh, poros) for tau, Dh, poros, a in zip(tau_cs, Dh_cs, poroses, expand_cores)]
    Aw_cores = [Aw_derive(alph,epsilon,d,N) for alph,epsilon,d,N in zip(alpha_cores,poroses,ds,Ns)]
    alpha_anns = [alpha_ann(cf1, cf2, rho, tau, mu, u, d, poros,D) for tau, d, poros in zip(tau_as, ds, poroses)]
    # alpha_anns = [alpha(cf1, cf2, rho, tau, mu, u, Dh, poros) for tau, Dh, poros, a in zip(tau_as, Dh_as, poroses, expand_cores)]
    Aw_anns = [Aw_derive(alph,epsilon,d,N) for alph,epsilon,d,N in zip(alpha_anns,poroses,ds,Ns)]

    Ergun_alphas = [Ergun_alpha(A_E, poros, d) for A_E, poros, d in zip(A_Es, poroses, ds)]

    beta_cores = [beta_new(cD, tau, l, poros) for tau, l, poros in zip(tau_cs, ls, poroses)]
    Bw_cores = [Bw_derive(bet,epsilon,d,N) for bet,epsilon,d,N in zip(beta_cores,poroses,ds,Ns)]
    beta_anns = [beta_new(cD, tau, l, poros) for tau, l, poros in zip(tau_as, ls, poroses)]
    Bw_anns = [Bw_derive(bet,epsilon,d,N) for bet,epsilon,d,N in zip(beta_anns,poroses,ds,Ns)]

    Ergun_betas = [Ergun_beta(B_E, poros, d) for B_E, poros, d in zip(B_Es, poroses, ds)]

    ann_fs = [2/N - 1/N**2 for N in Ns]
    weight_aws = [af*awa + (1-af)*awc for af,awa,awc in zip(ann_fs, Aw_anns, Aw_cores)]
    weight_bws = [af*bwa + (1-af)*bwc for af,bwa,bwc in zip(ann_fs, Bw_anns, Bw_cores)]


    fig1, ax1 = plt.subplots()
    Cheng_Aws = [A_E/M_factor(N,poros)**2 for A_E,N,poros in zip(A_Es,Ns,poroses)]
    ax1.plot(Ns,weight_aws,label="Flow paths (weight avg)")
    ax1.plot(Ns, Aw_cores, label="Core", linestyle='dashed')
    ax1.plot(Ns, Aw_anns, label="Annulus", linestyle='dashed')
    ax1.plot(Ns,Cheng_Aws,label="Cheng")
    ax1.set_xlabel("D/d ratio")
    ax1.set_ylabel("Aw")
    ax1.legend()
    plt.show()

    fig2, ax2 = plt.subplots()
    Cheng_Bws = [B_E / M_factor(N, poros) for B_E, N, poros in zip(B_Es, Ns, poroses)]
    ax2.plot(Ns,weight_bws,label="Flow paths (weight avg)")
    ax2.plot(Ns, Bw_cores, label="Core", linestyle='dashed')
    ax2.plot(Ns, Bw_anns, label="Annulus", linestyle='dashed')
    ax2.plot(Ns, Cheng_Bws,label="Cheng")
    ax2.set_xlabel("D/d ratio")
    ax2.set_ylabel("Bw")
    ax2.legend()
    plt.show()

    # fig3, ax3 = plt.subplots()
    # # ax3.plot(Ns, ls, label="layer heights")
    # # ax3.plot(Ns, ns, label="number per layer")
    # ax3.plot(Ns, tau_as, label="annular tortuosity ")
    # ax3.plot(Ns, tau_cs,label="core tortuosity")
    # ax3.set_xlabel("D/d ratio")
    # ax3.set_ylabel("quantity")
    # ax3.legend()
    # plt.show()


if __name__ == "__ain__":
    D = Length(1,Length.m)
    N = 2.7014
    d = D/N
    print("d:",d)
    n = 5
    n_a = na(N)
    print("na:",n_a)
    At = Atotal_core(D,d,n,n_a)
    print("AT:",At)
    print("max a:",np.pi*D**2/4)
    path_a = Apath_core(d)
    print("calcpath:",path_a)
    N_paths = At/path_a
    print(N_paths)
    N_guess = n-1
    Aguess = At/N_guess
    print("Aguess:",Aguess)
    a = np.sqrt((Aguess+np.pi*d**2/4)/d**2)
    print("increase factor:",a)

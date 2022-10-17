from dimension import *
import matplotlib.pyplot as plt


def tau_ann(l: Length, d: Length) -> float:
    """
    Tortuosity in the annulus. Calculated by averaging the path length when traveling directly upwards until meeting
    the sphere, then flowing along sphere for varying initial x positions.
    :param l: layer overlap (y difference between sphere centers in consecutive layers)
    :param d: sphere diameter
    :return: tau_a
    """
    return float(1 + d * (np.pi ** 2 - 8) / (8 * np.pi * l))


def tau_core(l: Length, d: Length) -> float:
    """
    Tortuosity in the core. Calculated by calculating the path length when traveling directly upwards until meeting
    the sphere, or starting flow along sphere and then flowing directly upwards.
    :param l: layer overlap (y difference between sphere centers in consecutive layers)
    :param d: sphere diameter
    :return: tau_c
    """
    return float(1 + d * (np.pi / 2 - 1) / l)


def A_intercept(D: Length, d: Length) -> Derived:
    """
    The area of the circle which intersects a larger circle
    :param d: sphere diameter (m)
    :param D: Inner diameter of the tube
    :return: A_int (m^2)
    """
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
    return d * np.sqrt((N ** 2 + 2 * N - n + na + 1 - 4 * A_int / (np.pi * d ** 2)) / (n-1))


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
            D_H ** (1 + cf2) * phi ** (1 - cf2))


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


def na(N: float) -> int:
    if N <= 2:
        return 1
    theta = np.arccos(1 / np.sqrt((N ** 2 - 2 * N + 1) / (N ** 2 - 2 * N)))
    return int(np.pi / theta)


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

# TODO: porosity calculation at core and annulus?

if __name__ == "__main__":
    D = Length(1, Length.inch)
    N_calcs = 100
    Ns = np.linspace(2.05, 2.95, N_calcs)
    ds = [D / N for N in Ns]
    ns_ls = [Guo_params(N, d) for N, d in zip(Ns, ds)]
    print(ns_ls)
    ns = [nl[0] for nl in ns_ls]
    print(ns)
    ls = [nl[1] for nl in ns_ls]
    [print(2 * l / d) for l, d in zip(ls, ds)]
    tau_as = [tau_ann(l, d) for l, d in zip(ls, ds)]
    tau_cs = [tau_core(l, d) for l, d in zip(ls, ds)]
    n_as = [na(N) for N in Ns]
    Dh_as = [Dh_ann(D, d, n_a) for d, n_a in zip(ds, n_as)]
    print("dh ann:",Dh_as)
    Dh_cs = [Dh_core(D, d, n, n_a) for d, n, n_a in zip(ds, ns, n_as)]
    print("dh core:", Dh_cs)
    poroses = [poros_Guo(N) for N in Ns]

    rho = Derived(998.2, units={Unit(Mass, Mass.kg): 1, Unit(Length, Length.m): -3})
    mu = Derived(0.001003,
                 units={Unit(Mass, Mass.kg): 1, Unit(Length, Length.m): -1, Unit(Time, Time.s): -1})

    vol = Volume(1, Volume.gal)
    duration = Time(1, Time.minute)
    mdot = rho * vol / duration
    area = np.pi * D ** 2 / 4
    u = mdot / (rho * area)

    cf1 = 64
    cf2 = 1
    cD = 0.62

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

    alpha_cores = [alpha(cf1, cf2, rho, tau, mu, u, Dh, poros) for tau, Dh, poros in zip(tau_cs, Dh_cs, poroses)]
    Aw_cores = [Aw_derive(alph,epsilon,d,N) for alph,epsilon,d,N in zip(alpha_cores,poroses,ds,Ns)]
    alpha_anns = [alpha(cf1, cf2, rho, tau, mu, u, Dh, poros) for tau, Dh, poros in zip(tau_as, Dh_as, poroses)]
    Aw_anns = [Aw_derive(alph,epsilon,d,N) for alph,epsilon,d,N in zip(alpha_anns,poroses,ds,Ns)]

    Ergun_alphas = [Ergun_alpha(A_E, poros, d) for A_E, poros, d in zip(A_Es, poroses, ds)]

    beta_cores = [beta(cD, tau, Dh, poros) for tau, Dh, poros in zip(tau_cs, Dh_cs, poroses)]
    Bw_cores = [Bw_derive(bet,epsilon,d,N) for bet,epsilon,d,N in zip(beta_cores,poroses,ds,Ns)]
    beta_anns = [beta(cD, tau, Dh, poros) for tau, Dh, poros in zip(tau_as, Dh_as, poroses)]
    Bw_anns = [Bw_derive(bet,epsilon,d,N) for bet,epsilon,d,N in zip(beta_anns,poroses,ds,Ns)]

    Ergun_betas = [Ergun_beta(B_E, poros, d) for B_E, poros, d in zip(B_Es, poroses, ds)]


    #
    # for i in range(N_calcs):
    #     print(Ns[i])
    #     # # print("alpha core:", alpha_cores[i])
    #     # # print("alpha ann:", alpha_anns[i])
    #     # print("alpha total:", alpha_anns[i]+alpha_cores[i])
    #     # print("alpha ergun:", Ergun_alphas[i])
    #     # # print("beta core:", beta_cores[i])
    #     # # print("beta ann:", beta_anns[i])
    #     # print("beta total:", beta_anns[i]+beta_cores[i])
    #     # print("beta ergun:", Ergun_betas[i])
    #     print("Aw me:",Aw_cores[i]+Aw_anns[i])
    #     print("Aw cheng:", A_Es[i]/M_factor(Ns[i],poroses[i])**2)
    #     print("Bw me:",Bw_cores[i]+Bw_anns[i])
    #     print("Bw cheng:", B_Es[i]/M_factor(Ns[i],poroses[i]))

    Cheng_Aws = [A_E/M_factor(N,poros)**2 for A_E,N,poros in zip(A_Es,Ns,poroses)]
    plt.plot(Ns,[sum(aws) for aws in zip(Aw_cores,Aw_anns)])
    plt.plot(Ns,Cheng_Aws)
    plt.show()

    Cheng_Bws = [B_E / M_factor(N, poros) for B_E, N, poros in zip(B_Es, Ns, poroses)]
    plt.plot(Ns,[sum(bws) for bws in zip(Bw_cores,Bw_anns)])
    plt.plot(Ns, Cheng_Bws)
    plt.show()
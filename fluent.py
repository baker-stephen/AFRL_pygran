import matplotlib.pyplot as plt
import numpy as np

from dimension import *

def poros_Foumeny(N):
    return .383 + .254*np.power(N,-0.923)*np.power(.723*N-1,-0.5)

def poros_Cheng(N):
    eps_1 = .8*(N-1)**.27
    eps_2 = .38*(1+(1/(N-1))**1.9)
    return (eps_1**(-3) + eps_2**(-3))**(-1/3)

def poros_Sato(N):
    #Used by Kerr thesis
    return 0.3494+0.4381/N

def poros_Guo(N):
    if N >= 3.0 or 1.866 < N < 2.0:
        # print("out of range for Guo")
        return -1

    k1=1
    k2=1

    if N<=1.866:
        k1 = N**(-2)
        k2 = 1/(np.sqrt(1 - (N-1)**2))
    elif 2<=N<=3:
        n = 2
        alpha = np.pi/4
        if 2.1547<=N<2.4142:
            n = 3
            alpha = np.pi/6
        elif 2.4142<=N<2.7013:
            n = 4
            alpha = 22.5 * (np.pi/180)
        elif 2.7013<=N:
            n = 5
            alpha = 18 * (np.pi/180)
        k1 = n/(N**2)
        #Guo says there should be a Dp in there, but units don't work??
        k2 = 1/np.sqrt(1-((N-1)*np.sin(alpha))**2)

    return 1-(2/3)*k1*k2

def A_E_B_E(porosity,N,a1,a2,b1,b2,m):
    A_E = (a1 + (a2*porosity/(1-porosity))*(N/(N-1))**m)
    B_E = (b1*np.power((1-porosity)/porosity,1/3)+b2*(N/(N-1))**m)
    return A_E,B_E

def AwBw(porosity,N, A_E, B_E):
    M = 1 + 2/(3*N*(1-porosity))
    print("M: ",M)
    Aw = A_E / M**2
    Bw = B_E / M
    return Aw,Bw

def C1C2_old(A_E, B_E, porosity, Dp):
    C1 = A_E * (1 - porosity) ** 2 / (porosity ** 3 * Dp ** 2)
    C2 = 2 * B_E * (1 - porosity) / (porosity ** 3 * Dp)
    return C1, C2

def C1C2(A_E, B_E, porosity, Dp):
    C1 = A_E * (1 - porosity) ** 2 / (porosity ** 2 * Dp ** 2)
    C2 = 2 * B_E * (1 - porosity) / (porosity * Dp)
    return C1, C2

def find_mdot(density: Derived, h: Length, C_D: float, diameter: Length) -> Derived:
    LA_rad = Length(6371.471,Length.km)
    lil_g = Fun.m_Earth*Fun.G/LA_rad**2
    vel = C_D * np.sqrt(2 * lil_g * h)
    return density*vel*np.pi*(diameter/2)**2

def cheng_dp(A_E:float, B_E:float, porosity:float, density:Derived, dyn_viscosity:Derived, part_diameter:Length, superficial_velo:Derived) -> Derived:
    return A_E * (1 - porosity) ** 2 * dyn_viscosity * superficial_velo / (porosity ** 3 * part_diameter ** 2) + B_E * (1 - porosity) * density * superficial_velo ** 2 / (porosity ** 3 * part_diameter)

def reynold(porosity:float, part_diameter:Length,superficial_velo:Derived,kinematic_visosity:Derived) -> float:
    return float(superficial_velo*part_diameter/(kinematic_visosity*(1-porosity)))

def turb_intensity(reynold:float) -> float:
    """
    Ansys guide 7.4.2.1.3, ref 165
    :param reynold: Reynolds number (dimless)
    :return: turbulent intensity I
    """
    return .16*reynold**(-1/8)


if __name__ == "__main__":
    #Constants
    a1 = 185
    a2 = 17
    b1 = 1.3
    b2 = .03
    m = 2

    rho = Derived(998.2,
                  numerator=[Unit(Mass,Mass.kg)],
                  denominator={Unit(Length,Length.m):3})
    dyn_viscosity = Derived(0.001003,
                            numerator=[Unit(Mass,Mass.kg)],
                            denominator=[Unit(Length,Length.m),Unit(Time,Time.s)])
    water_height = Length(7,Length.inch)
    discharge_coef = .62
    leng = Length(12,Length.inch)

    #Parameters to vary
    Dp_mini = Length(1,Length.mm)
    Dp_super_small = Length(1/8,Length.inch)
    Dp_small = Length(1/4,Length.inch)
    Dp_mid = Length(3/8,Length.inch)
    Dp_large = Length(7/16,Length.inch)
    D_one_quarter = Length(.26,Length.inch)
    D_one_quarter_outside = Length(.375, Length.inch)
    D_half = Length(5 / 8, Length.inch)
    D_half_outside = Length(13 / 16, Length.inch)
    D_three_quarter = Length(7 / 8, Length.inch)
    D_three_quarter_outside = Length(17 / 16, Length.inch)
    D_one = Length(1.029, Length.inch)
    D_one_outside = Length(1.315, Length.inch)
    D_one_one_half = Length(1.59, Length.inch)
    D_one_one_half_outside = Length(1.9, Length.inch)

    #Set this experiment parameters here
    Dp = Length(1/4, Length.inch)
    D = Length(0.602, Length.inch)

    print("water mass flow rates from experimental data")
    rho_h20_exp = Derived(998.23, numerator=[Unit(Mass,Mass.kg)],
                          denominator={Unit(Length,Length.m):3})
    vols = [Volume(v,Volume.gal) for v in np.linspace(0.01,3.00,300)]
    vol = Volume(1.88,Volume.gal)
    duration = Time(60,Time.s)
    mass_flow_exp = vol*rho_h20_exp/duration
    mdots_list = [v*rho_h20_exp/duration for v in vols]
    print("mass_flow_exp:",mass_flow_exp)
    print(float(mass_flow_exp))
    print("")

    #outputs
    print("porous zone inputs:")
    N = float(D/Dp)
    print("N:",N)
    porosity = poros_Foumeny(N)
    # porosity = poros_Cheng(N)
    print("Foumeny poros:", porosity)
    porosity = 0.466470357313992 #TODO Manually enter Pygran fitted average porosity
    # print("pygran porosity: ",porosity)
    print("cheng poros:",poros_Cheng(N))
    print("sato poros:",poros_Sato(N))
    vol_to_fill = (np.pi * (D/2)**2 * leng) * (1-porosity)
    vol_sphere = (4/3) * np.pi * (Dp/2) ** 3
    N_spheres = vol_to_fill / vol_sphere
    print("N spheres:",N_spheres)
    # porosity = 0.4271345852237037
    #pygran results:
    # num_spheres = 73
    # V_spheres = num_spheres*(4/3)*np.pi*(Dp/2)**3
    # pygran_len = Length(4,Length.inch)
    # V_cyl = np.pi*(D/2)**2*pygran_len
    # porosity = (V_cyl-V_spheres)/V_cyl

    A_E, B_E = A_E_B_E(porosity,N,a1,a2,b1,b2,m)
    print("A_E: ",A_E)
    print("B_E: ",B_E)
    Aw, Bw = AwBw(porosity,N,A_E,B_E)
    print("Aw: ",Aw)
    print("Bw: ",Bw)
    C1, C2 = C1C2(A_E, B_E, porosity, Dp)
    print("viscous resistance (C1): ",C1)
    print(float(C1))
    print("inertial resistance (C2): ",C2)
    print(float(C2))
    C1_old, C2_old = C1C2_old(A_E, B_E, porosity, Dp)
    print("C1, C2 old formulation")
    print("C1 old: ", C1_old)
    print("C2 old: ", C2_old)
    print("")

    print("Inlet/outlet conditions:")
    # mdot = find_mdot(rho,water_height,discharge_coef,D)
    mdot = mass_flow_exp
    print(mdot.get_special([Unit(Volume, Volume.gal)], warn=False))
    print("mdot: ",mdot)
    print(float(mdot))
    area = np.pi*(D/2)**2
    superficial_velo = mdot / (rho * area)
    print("OG U:",superficial_velo)
    # superficial_velo = Derived(4.2725326,numerator=[Unit(Length,Length.m)],
    #                       denominator=[Unit(Time,Time.s)])*porosity #TODO: manually entered velo from fluent
    poroses = np.linspace(0.466470357313992*.98,0.466470357313992*1.02,500)
    dpdxes = [cheng_dp(A_E,B_E,p,rho,dyn_viscosity,Dp,superficial_velo) for p in poroses]
    plt.plot(poroses,dpdxes)
    plt.show()
    dps = []
    for md in mdots_list:
        U = md / (rho * area)
        dps.append(cheng_dp(A_E, B_E, porosity, rho, dyn_viscosity, Dp, U) * leng)
    # fig, ax = plt.subplots()
    # ax.set_xlabel('Volumetric Flow Rate (gal/min)')
    # ax.set_ylabel('Total Pressure Drop (psid)')
    # ax.plot([md.get_special([Unit(Volume,Volume.gal)], warn=False)[0]*0.06 for md in mdots_list], [dp.get_special(units=[Unit(Pressure,Pressure.psi)], warn=False)[0] for dp in dps])
    # plt.show()
    first_cheng = cheng_dp(A_E,B_E, porosity, rho, dyn_viscosity, Dp, superficial_velo)
    print("First cheng: ",first_cheng)
    print(first_cheng*Length(1,Length.inch))
    print((first_cheng*Length(1,Length.inch)).get_special([Unit(Pressure,Pressure.psi)]))
    re = reynold(porosity,Dp,superficial_velo,dyn_viscosity/rho)
    print("reynold: ",re)
    turb = turb_intensity(re)
    print("turb %: ",turb*100)
    print("D_H: ",D)
    print(float(D))


    #checking fluent derived flow velocity:
    print("")
    print("Fluent output checking:")
    v_fluent = Derived(1.161137,numerator=[Unit(Length,Length.m)],denominator=[Unit(Time,Time.s)])
    r_fluent = np.sqrt(mdot / (rho * np.pi * v_fluent))
    print("fluent: ",r_fluent.get_special([Unit(Length,Length.inch)]))
    # For straight inlet testing:
    print("original: ",(D/2).get_special([Unit(Length,Length.inch)]))

    #calculate pressure drop via Cheng:
    print("U: ",superficial_velo)
    print(float(superficial_velo))
    print("physical avg v: ",superficial_velo/porosity)
    dpdx = cheng_dp(A_E,B_E,porosity,rho,dyn_viscosity,Dp,superficial_velo)
    print("dpdx: ",dpdx)
    print("dpdx: ",dpdx.get_special([Unit(Pressure,Pressure.psi)]))
    print("dpdx full:",float(dpdx))
    dp_total = dpdx*leng
    print("dp_total: ", dp_total.get_special(units=[Unit(Pressure,Pressure.psi)]))
    print(float(dp_total))

    # 3/4" D, 1/4" Dp, straight pipe only -
    # dp_fluent_laminar = Pressure(3.604428888E+05 - 9.717660602E+00,Pressure.Pa) #For fine mesh, laminar case
    # dp_fluent_turbulent = Pressure(3.607939438E+05 - 2.801655739E+01,Pressure.Pa) #For fine mesh, turbulent
    # print("dp laminar: ",dp_fluent_laminar)
    # print("dp turbulent: ", dp_fluent_turbulent)
    # print("percent diff: ",abs(dp_fluent_laminar-dp_fluent_turbulent)/dp_fluent_turbulent)







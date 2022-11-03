from dimension import *

def poros_Foumeny(N):
    #Great accuracy at slightly higher N
    return .383 + .254*np.power(N,-0.923)*np.power(.723*N-1,-0.5)

def poros_Cheng(N):
    #Solid accuracy for the widest range of N
    eps_1 = .8*(N-1)**.27
    eps_2 = .38*(1+(1/(N-1))**1.9)
    return (eps_1**(-3) + eps_2**(-3))**(-1/3)

def poros_Sato(N):
    #Used by Kerr thesis
    return 0.3494+0.4381/N

def poros_Guo(N):
    #Analytically derived. Should be exact for perfectly settled/ordered spheres. Very limited range of applicability.
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
    # Generally unnecessary. Use A_E, B_E
    M = 1 + 2/(3*N*(1-porosity))
    print("M: ",M)
    Aw = A_E / M**2
    Bw = B_E / M
    return Aw,Bw

def C1C2_old(A_E, B_E, porosity, Dp):
    # Assuming superficial velocity = physical
    C1 = A_E * (1 - porosity) ** 2 / (porosity ** 3 * Dp ** 2)
    C2 = 2 * B_E * (1 - porosity) / (porosity ** 3 * Dp)
    return C1, C2

def C1C2(A_E, B_E, porosity, Dp):
    # Converting superficial (U) to physical (v) by v = U/porosity.
    C1 = A_E * (1 - porosity) ** 2 / (porosity ** 2 * Dp ** 2)
    C2 = 2 * B_E * (1 - porosity) / (porosity * Dp)
    return C1, C2

def find_mdot(density: Derived, h: Length, C_D: float, diameter: Length) -> Derived:
    # No longer used. This was for when we were using buckets as source.
    # This function determined the flow rate out of the bucket.
    LA_rad = Length(6371.471,Length.km)
    lil_g = Fun.m_Earth*Fun.G/LA_rad**2
    vel = C_D * np.sqrt(2 * lil_g * h)
    return density*vel*np.pi*(diameter/2)**2

def cheng_dp(A_E:float, B_E:float, porosity:float, density:Derived, dyn_viscosity:Derived, part_diameter:Length, superficial_velo:Derived) -> Derived:
    # This is the big one. dP/dx as given by Cheng.
    return A_E * (1 - porosity) ** 2 * dyn_viscosity * superficial_velo / (porosity ** 3 * part_diameter ** 2) + B_E * (1 - porosity) * density * superficial_velo ** 2 / (porosity ** 3 * part_diameter)

def reynold(porosity:float, part_diameter:Length,superficial_velo:Derived,kinematic_visosity:Derived) -> float:
    # Reynolds number in porous media.
    return float(superficial_velo*part_diameter/(kinematic_visosity*(1-porosity)))

def turb_intensity(reynold:float) -> float:
    """
    Ansys guide 7.4.2.1.3, ref 165
    :param reynold: Reynolds number
    :return: turbulent intensity I
    """
    return .16*reynold**(-1/8)


if __name__ == "__main__":
    #Empirically derived Cheng constants
    a1 = 185
    a2 = 17
    b1 = 1.3
    b2 = .03
    m = 2

    # Density of water defined in Fluent materials.
    rho = Derived(998.2, {Unit(Mass, Mass.kg): 1, Unit(Length, Length.m): -3, })

    # Dynamic viscosity of water defined in Fluent materials.
    dyn_viscosity = Derived(0.001003,{Unit(Mass, Mass.kg): 1, Unit(Length, Length.m): -1,Unit(Time,Time.s):-1})

    # Pipe length
    leng = Length(12,Length.inch)

    #Parameters to vary
    Dp_mini = Length(1,Length.mm)
    Dp_super_small = Length(1/8,Length.inch)
    Dp_small = Length(1/4,Length.inch)
    Dp_mid = Length(3/8,Length.inch)
    Dp_large = Length(7/16,Length.inch)
    D_one_quarter = Length(.26,Length.inch)
    D_one_quarter_outside = Length(.375, Length.inch)
    D_half = Length(0.602, Length.inch)
    D_half_outside = Length(13 / 16, Length.inch)
    D_three_quarter = Length(7 / 8, Length.inch)
    D_three_quarter_outside = Length(17 / 16, Length.inch)
    D_one = Length(1.029, Length.inch)
    D_one_outside = Length(1.315, Length.inch)
    D_one_one_half = Length(1.59, Length.inch)
    D_one_one_half_outside = Length(1.9, Length.inch)

    # TODO: Set the sphere (pellet) diameter
    Dp = Length(2.8,Length.mm)
    # TODO: Set he inner diameter
    D = D_half
    N = float(D / Dp)
    print("N:",N)
    empirical_porosity = poros_Foumeny(N) if poros_Guo(N) == -1 else poros_Guo(N)
    # empirical_porosity = poros_Cheng(N)
    # TODO: Manually enter Pygran fitted volume-averaged porosity from outputs.txt
    # porosity = 0.416128576115776
    porosity = empirical_porosity
    D_H = np.sqrt(8*Dp**2*porosity**3/(9*(1-porosity)**2))
    print("D_H:",D_H.get_special([Unit(Length,Length.inch)]))
    A_H_ratio = D_H**2/D**2
    print("A_H_ratio:",A_H_ratio)

    # Calulate mass flow rate from volumetric flow rate of experiment
    rho_h20_exp = Derived(998.23, {Unit(Mass, Mass.kg): 1, Unit(Length, Length.m): -3, })

    # TODO: adjust gallons per minute if necessary/desired
    vol = Volume(1,Volume.gal)
    duration = Time(1,Time.minute)
    mass_flow_exp = vol*rho_h20_exp/duration
    mdot = mass_flow_exp

    #outputs

    print("Mass flow inlet/outlet boundary conditions:")
    print("mass flow rate:",mdot)
    print(float(mdot))
    area = np.pi*(D/2)**2
    superficial_velo = mdot / (rho * area)
    Re = reynold(porosity, Dp, superficial_velo, dyn_viscosity / rho)
    print("Turbulent Intensity:",turb_intensity(Re)*100)
    print("Hydraulic Diameter:",float(D))

    print("\nNamed Expressions:")
    print("MAKE SURE THESE ARE UPDATED BETWEEN RUNS")
    print("Dp:", float(Dp))
    print("ID:", float(D))
    print("poros_fourier:", "COPY FROM outputs.txt, the line under 'Fourier fit degree __ -'")

    print("\nThe rest can be safely copy/pasted between runs:")
    print("a1:", a1)
    print("a2:", a2)
    print("b1:", b1)
    print("b2:", b2)
    print("m:", m)
    print("N:", "ID / Dp")
    print("p_r:", "sqrt(x^2+y^2)")
    print("poros:", "max(min(poros_fourier,1),0)")
    print("A_E:", "(a1 + (a2*poros/(1-poros))*(N/(N-1))^2)")
    print("B_E:", "(b1*((1-poros)/poros)^(1/3)+b2*(N/(N-1))^2)")
    # Adjusted C1, C2. Use commented lines for original.
    # Did some testing
    # print("C1:", "A_E * (1 - poros) ^ 2 / (poros ^ 2 * Dp ^ 2)")
    print("C1:", "A_E * (1 - poros) ^ 2 / (poros ^ 3 * Dp ^ 2)") # One extra power of porosity in denominator
    # print("C2:", "2*B_E*(1-poros)/(poros*Dp)")
    print("C2:", "2*B_E*(1-poros)/(poros**3*Dp)") # Two extra powers of porosity in denominator

    print("\n1D estimates for comparison:")
    print("physical velocity:",superficial_velo/porosity)
    A_E, B_E = A_E_B_E(porosity, N, a1, a2, b1, b2, m)
    print("A_E: ",A_E)
    print("B_E: ",B_E)
    C1, C2 = C1C2_old(A_E, B_E, porosity, Dp)
    print("viscous resistance (C1): ",C1)
    print(float(C1))
    print("inertial resistance (C2): ",C2)
    print(float(C2))
    dpdx = cheng_dp(A_E,B_E,porosity,rho,dyn_viscosity,Dp,superficial_velo)
    print("dpdx: ",dpdx)
    print("dpdx: ",dpdx.get_special([Unit(Pressure,Pressure.psi)]))
    print("dpdx (Pa/m):", float(dpdx))
    dp_total = dpdx*leng
    print("dp_total: ", dp_total.get_special([Unit(Pressure,Pressure.psi)]))
    print("dp_total (Pa): ",float(dp_total))

    print("\nother possibly relevant outputs:")
    print("N:", N)
    print("Foumeny poros:", poros_Foumeny(N))
    print("Cheng poros:", poros_Cheng(N))
    print("Guo poros:", poros_Guo(N))
    print("Sato poros:",poros_Sato(N))
    C1_other, C2_other = C1C2_old(A_E, B_E, porosity, Dp)
    print("C1 other: ", C1_other)
    print(float(C1_other))
    print("C2 other: ", C2_other)
    print(float(C2_other))

    sdpdx = Derived(20.523858267716538,{Unit(Pressure,Pressure.Pa):1,Unit(Length,Length.m):-1})
    RHS = sdpdx/dyn_viscosity
    print("RHS:",RHS)
    print(float(RHS))
    mat_Q = mdot/rho
    print("mat_Q:",mat_Q)
    tau = 1.54361554933
    dhx = Dp*0.408248290464
    Recivan = rho*tau*superficial_velo*dhx/(porosity*dyn_viscosity)
    print("REcivan:",Recivan)






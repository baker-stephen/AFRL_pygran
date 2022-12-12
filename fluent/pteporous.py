from dimension import *

#Homework 1

def kozeny_carmen(Gamma: Length, phi: float)->Permeability:
    """
    PMTP Eqn 2-7
    :param Gamma: Intrinsic parameter
    :param phi: porosity
    :return: K (permeability, m^2)
    """
    return Permeability(float(Gamma**2*phi**2/(1-phi)**3),Permeability.m2)

def VTF(E: TempEnergy, T: Temperature, T_c: Temperature, T_0: Temperature, K_0: Permeability)->Permeability:
    """
    PMTP p66 eqn 2.150
    :param E: activation energy
    :param T: temperature
    :param T_c: characteristic limit temperature
    :param T_0: reference temperature
    :param K_0: reference permeability
    :return: modified permeability
    """
    return K_0*np.exp((-E/Fun.R)*(1/(T-T_c) - 1/(T_0-T_c)))

def grain_shape(V: Volume, A: Derived)->float:
    """
    RFD ch3 p11 eqn3.27
    :param V: grain volume
    :param A: grain surface area
    :return: sphericity
    """
    return float(6**(2/3)*np.pi**(1/3)*V**(2/3)/A)

#Midterm 2

def Klink(K_inf: Permeability, b: Pressure, p: Pressure)->Permeability:
    return Permeability(float(K_inf*(1+b/p)),Permeability.m2)

def Darcy_velo(K: Permeability, mu: Viscosity, dPdx: Derived)->Derived:
    return -K*dPdx/mu

def Darcy_drop(K: Permeability, mu: Viscosity, u: Derived)->Derived:
    return -u*mu/K

if __name__ == "__main__":

    # sig = Derived(1,{Unit(Force,Force.dyn):1,Unit(Length,Length.cm):-1})
    # rc = Length(1,Length.mm)
    # Pc = Pressure(1,Pressure.psi)
    # out = 2*sig/(rc*Pc)
    # print(out)

    z1 = 0.86
    z2 = 0.14
    k1 = 1.2
    k2 = 0.38
    v_1 = -(z1/(k2-1)+z2/(k1-1))
    print(v_1)
    a = (k1-1)*(k2-1)
    b = -z1*(k2-1)-z2*(k1-1)+k1+k2-2
    c = 1-z1-z2
    v_2a = (-b+np.sqrt(b**2-4*a*c))/(2*a)
    print(v_2a)
    v_2b = (-b-np.sqrt(b**2-4*a*c))/(2*a)
    print(v_2b)
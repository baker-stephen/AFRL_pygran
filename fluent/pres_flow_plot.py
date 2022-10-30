import numpy as np

from dimension import *
import matplotlib.pyplot as plt
from fluent import poros_Foumeny, poros_Guo, cheng_dp, A_E_B_E

if __name__ == "__main__":
    a1 = 185
    a2 = 17
    b1 = 1.3
    b2 = .03
    m = 2

    Dp = Length(2.8,Length.mm)
    D = Length(0.602, Length.inch)
    area = np.pi * (D / 2) ** 2
    dyn_viscosity = Derived(0.001003,{Unit(Mass,Mass.kg):1,Unit(Length,Length.m):-1,Unit(Time,Time.s):-1})
    density = Derived(998.23,{Unit(Mass,Mass.kg):1,Unit(Length,Length.m):-3})
    N = float(D/Dp)
    porosity_base = poros_Foumeny(N) if poros_Guo(N) == -1 else poros_Guo(N)
    fig, ax = plt.subplots()
    Vdots = [0.36, .5, 1.035, 1.14, 1.5725, 1.92, 2.10585, 2.63]
    with open('2pt8_25pt4-0pt602_dpdx.csv','w') as f:
        f.write(' ')
        for vdot in Vdots:
            f.write(', '+str(vdot))
        f.write('\n')
        for factor in [0.9,1,1.1]:
            porosity = porosity_base*factor
            f.write(str(porosity))
            A_E, B_E = A_E_B_E(porosity, N, a1, a2, b1, b2, m)
            # numpoints = 100
            numpoints = len(Vdots)
            presses = np.zeros(numpoints)
            # Vdots = np.linspace(0.3, 2.7, numpoints)
            for i in range(numpoints):
                Vdot = Derived(Vdots[i], {Unit(Volume, Volume.gal):1, Unit(Time, Time.minute):-1})
                superficial_velo = Vdot / area
                # print(superficial_velo)
                presses[i] = cheng_dp(A_E,B_E,porosity,density,dyn_viscosity,Dp,superficial_velo)
                f.write(', '+str(presses[i]))

            f.write('\n')
            ax.scatter(Vdots, presses,label="Por="+str(porosity))

    ax.legend()
    ax.set_ylabel("Pressure Drop (Pa/m)")
    ax.set_xlabel("Volumetric Flow Rate (gal/min)")
    plt.show()
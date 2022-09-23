from dimension import *
from fluent import poros_Foumeny, poros_Guo

if __name__=="__main__":
    D = Length(0.602, Length.inch)
    Dp = Length(1/4, Length.inch)
    N = D/Dp
    leng = Length(12,Length.inch)
    porosity = poros_Foumeny(N) if poros_Guo(N) == -1 else poros_Guo(N)
    vol_to_fill = (np.pi * (D/2)**2 * leng) * (1-porosity)
    vol_sphere = (4/3) * np.pi * (Dp/2) ** 3
    N_spheres = vol_to_fill / vol_sphere
    print("N spheres:",N_spheres)
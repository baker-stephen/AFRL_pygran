import numpy as np

# Define constants

e_charge = 1.6022e-19  # C

R_gas = 8.3145e3  # J / K * kmol
# reference temperature
T_f = 298.15  # K

R_E = 6378.14e3  # m
R_E_polar = 6356.8e3  # m
R_E_V = 6371.0e3  # m
R_GEO = 42164.2e3  # m

G = 6.672e-11  # m^3 / kg * s^2
g_E_0 = 9.80665  # m/s^2

mu_E = 3.986e14  # m^3/s^2
M_E = 5.974e24  # kg
J_2_E = .00108263

mu_sun = 1.3271e20  # m^3/s^2
M_sun = M_E * 3.33e5

mu_Mercury = 2.18e13  # m^3/s^2
M_Mercury = .053 * M_E

mu_Venus = 3.249e14  # m^3/s^2
M_Venus = .815 * M_E
J_2_Venus = .000027

mu_Mars = 4.293e13
M_Mars = M_E * .107
J_2_Mars = .001964

F_s = 1363  # W/m^2

# avocado's
N_A = 6.022e26  # 1/kmol
# Boltzman
k_B = 1.3807e-23  # J / K

# electron mass
m_e = 9.1093837015e-31

# proton mass
m_p = 1.67262192369e-27

# speed of light
c = 2.9979e8  # m/s

# Stefan-Boltzmann constant
sigma_SB = 5.6703e-8  # W/m^2*K^4

# Sound pressure level reference
P_REF = 2e-5  # Pa

# permitivity of free space
epsilon_0 = 8.85e-12  # 8.85418782e-12 C^2/N*m^2

#Planck's constant
h_planck = 6.6261e-34 #J/K

#obliquity of the ecliptic
epsilon_obliquity = (23 + 26/60 + 22/3600) * np.pi / 180.0 #radians

#Sidereal day
day_sidereal = 86164.09 #seconds
year_sidereal = 31558149.504 #seconds
#from equinox to equinox
year_tropical = 31556925.19008 #s
#perihelion to perihelion
year_anomalistic = 31558432.896 #s

eclipses_day_LEO = 15
eclipses_year_GEO = 90

# Conversions
def inch_cm(inches):
    return inches * 2.54


def cm_m(cms):
    return cms / 100.0


def inch_m(inches):
    return cm_m(inch_cm(inches))


def mil_inch(mils):
    return mils * np.power(10.0, -3)


def mil_m(mils):
    return inch_m(mil_inch(mils))


def ft_m(feet):
    return feet * .3048


def mile_m(miles):
    return miles * 1609.344


def nmile_m(nmiles):
    return nmiles * 1852.0


def oz_kg(ounces):
    return ounces * .02835


def lb_kg(lbs):
    return lbs * 0.4536


def lbf_N(lbfs):
    return lbfs * 4.448


def slug_kg(slugs):
    return slugs * 14.59

def dyn_N(dyns):
    return dyns*1e-5


def bar_Pa(bars):
    return bars * np.power(10.0, 5)


def atm_Pa(atms):
    return atms * 101325.0


def psi_Pa(psis):
    return psis * 6894.76


def torr_Pa(torrs):
    return torrs * 133.3223684211

def liter_m3(liters):
    return liters * 1e-3

def gal_m3(gals):
    return gals * 0.00378541


def calth_J(cals):
    return cals*4.184

def kcal_J(kcals):
    return 4184*kcals

def eV_J(eVs):
    return eVs * e_charge


def eV_K(eVs):
    return eV_J(eVs) / k_B


def Au_m(Aus):
    return Aus * 1.496e11


def pc_m(pcs):
    return pcs * 3.09e16


def ly_m(lys):
    return lys * 9.46e15


def ydhms_sec(years=0, days=0, hours=0, mins=0, secs=0):
    """
    :return: total seconds
    :rtype: float
    """
    days = days + years * 365.242
    hours = hours + days * 24
    mins = mins + hours * 60
    secs = secs + mins * 60
    return secs


def deg_rad(degs):
    return degs * np.pi / 180.0


def rad_deg(rads):
    return rads * 180.0 / np.pi

def dms_rad(deg = 0.0, amin = 0.0, asec = 0.0):
    amin += asec/60
    deg += amin/60
    return deg_rad(deg)


def h_r(h):
    return h + R_E

def R_K(rankines):
    return rankines*1.8

def F_K(fahrens):
    return (fahrens-32)*(5/9) + 273.15

def F_R(fahrens):
    return fahrens + 459.67

def C_K(celsius):
    return celsius + 273.15


def K_C(kelvins):
    return kelvins - 273.15


def amu_kg(amus):
    return amus * 1.6605e-27


def mega_(megas):
    return 10 ** 6 * megas


def giga_(gigas):
    return 10 ** 9 * gigas


def milli_(millis):
    return 10 ** (-3) * millis


def kilo_(kilos):
    return 10 ** 3 * kilos

def micro_(micros):
    return 10**(-6) * micros

def nano_(nanos):
    return 10**(-9) * nanos

def dB_gain(dBs):
    return np.power(10.0, dBs / 10.0)


def gain_dB(gains):
    return 10 * np.log10(gains)

def Gauss_Tesla(gausses):
    return gausses*1e-4

def gamma_Tesla(gammas):
    return gammas*1e-9

def sfu_SI(sfus):
    """
    :param sfus: solar flux units
    :return: W/ m^2*Hz
    """
    return sfus*1e-22

def ergs_J(ergs):
    return ergs*1e-7

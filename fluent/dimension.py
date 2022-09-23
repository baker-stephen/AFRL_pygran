import constants_conversions as cc
from abc import ABC, abstractmethod
from collections import Counter
import numpy as np

#TODO: convert into derived denominator
#TODO: make all fundamental constants dimensional
#TODO: frequency (Hz)
#TODO: amount of substance (mol)
#TODO: Luminous intensity (cd)
#TODO: test charge
#TODO: test angles
#TODO: test pressure
#TODO: test force
#TODO: test energy
#TODO: test standalone Temperature
#TODO: test Time
#TODO: non-int power?
#TODO: other comparisons
#TODO: apply is comparable method to add, sub, etc.


class Dimension(ABC):
    @abstractmethod
    def get_default(self):
        pass

    def is_comparable(self, other) -> bool:
        if not other.__class__.mro().__contains__(Dimension):
            if (other.__class__ is float or other.__class__ is int or other.__class__ is np.float64) and \
                    self.is_dimless():
                return True
            raise Exception("Must both be of type dimensional or have no dimension to compare.")
        numCpy = self.numerator.copy()
        denCpy = self.denominator.copy()
        otherNum = other.numerator.copy()
        otherDenom = other.denominator.copy()
        if Counter(numCpy) != Counter(otherNum) or Counter(denCpy) != Counter(otherDenom):
            raise Exception("Both values must have same dimension to compare.")
        return True

    def __str__(self):
        ret_str = "{:e} "
        ret_str += self.unit_to_str(self.numerator,self.denominator)
        return ret_str.format(self.value)

    def __repr__(self):
        return self.__str__()

    def unit_to_str(self, numerator: list, denominator: list):
        ret_str = ""
        numCtr = Counter(numerator)
        denCtr = Counter(denominator)

        for key in numCtr:
            if numCtr[key] > 1:
                ret_str += key.minor+"^"+str(numCtr[key])+"*"
            else:
                ret_str += key.minor + "*"
        if len(self.numerator) > 0:
            ret_str = ret_str[0:len(ret_str) - 1]
        if len(self.denominator) > 0:
            ret_str += "/"
            for key in denCtr:
                if denCtr[key] > 1:
                    ret_str += key.minor+"^"+str(denCtr[key])+"*"
                else:
                    ret_str += key.minor + "*"
            ret_str = ret_str[0:len(ret_str) - 1]
        return ret_str

    def __mul__(self, other) -> 'Derived':
        newNum = self.numerator.copy()
        newDenom = self.denominator.copy()

        if not other.__class__.mro().__contains__(Dimension):
            if other.__class__ is int or other.__class__ is float or other.__class__ is np.float64:
                other = float(other)
                return Derived(self.value*other, numerator=newNum, denominator=newDenom)
            else:
                print("trying to multiply other: ",other.__class__)
                raise TypeError("Can only multiply by float, int, and other Dimension.")

        otherNumCopy = other.numerator.copy()
        otherDenomCopy = other.denominator.copy()

        i = 0
        while i < len(newNum):
            if otherDenomCopy.__contains__(newNum[i]):
                otherDenomCopy.remove(newNum[i])
                newNum.remove(newNum[i])
            else:
                i+=1

        i = 0
        while i < len(newDenom):
            if otherNumCopy.__contains__(newDenom[i]):
                otherNumCopy.remove(newDenom[i])
                newDenom.remove(newDenom[i])
            else:
                i += 1

        newNum.extend(otherNumCopy)
        newDenom.extend(otherDenomCopy)
        return Derived(self.value*other.value, numerator=newNum, denominator=newDenom)

    def __float__(self) -> float:
        #TODO: check if dimensional and warn
        return float(self.value)

    def __add__(self, other) -> 'Derived':
        if not other.__class__.mro().__contains__(Dimension):
            if (other.__class__ is float or other.__class__ is int) and \
                    self.is_dimless():
                return Derived(self.value+other)
            raise Exception("Must both be of type dimensional or have no dimension to add.")
        numCpy = self.numerator.copy()
        denCpy = self.denominator.copy()
        otherNum = other.numerator.copy()
        otherDenom = other.denominator.copy()
        if Counter(numCpy) != Counter(otherNum) or Counter(denCpy) != Counter(otherDenom):
            raise Exception("Both values must have same dimension to add.")

        return Derived(self.value + other.value,numerator=numCpy,denominator=denCpy)

    def __pow__(self, power, modulo=None) -> 'Derived':

        if modulo is not None:
            raise ValueError("Modulo of dimensional quantities not supported.")


        numCtr = Counter(self.numerator)
        denCtr = Counter(self.denominator)

        newNum = []
        newDen = []

        for key in numCtr:
            count = numCtr[key]
            count *= power
            if abs(count - int(count)) > 1e-5:
                raise ValueError("count for " + str(key) + " is " + str(count)+"\n"
                                 + "Non-integer powers of dimensions not supported.")
            count = int(count)
            for i in range(count):
                newNum.append(key)

        for key in denCtr:
            count = denCtr[key]
            count *= power
            if abs(count - int(count)) > 1e-5:
                raise ValueError("count for " + str(key) + " is " + str(count) + "\n"
                                 + "Non-integer powers of dimensions not supported.")
            count = int(count)
            for i in range(count):
                newDen.append(key)

        return Derived(pow(self.value, power), numerator=newNum, denominator=newDen)

    def __truediv__(self, other) -> 'Derived':
        newNum = self.numerator.copy()
        newDenom = self.denominator.copy()

        if not other.__class__.mro().__contains__(Dimension):
            if other.__class__ is int or other.__class__ is float or other.__class__ is np.float64:
                other = float(other)
                return Derived(self.value / other, numerator=newNum, denominator=newDenom)
            else:
                raise TypeError("Can only divide float, int, other Dimension.")

        otherNumCopy = other.numerator.copy()
        otherDenomCopy = other.denominator.copy()
        inverse = Derived(1/other.value, numerator=otherDenomCopy, denominator=otherNumCopy)
        return self * inverse

    def __sub__(self, other) -> 'Derived':
        if not other.__class__.mro().__contains__(Dimension):
            if (other.__class__ is float or other.__class__ is int) and \
                    (self.numerator is None or len(self.numerator) == 0) and \
                    (self.denominator is None or len(self.denominator) == 0):
                return Derived(self.value + other)
            raise Exception("Must both be of type dimensional or have no dimension to subtract.")

        return self + (-1*other)

    def __rmul__(self, other) -> 'Derived':
        return self.__mul__(other)

    def __rtruediv__(self, other) -> 'Derived':

        newNum = self.numerator.copy()
        newDenom = self.denominator.copy()

        if not other.__class__.mro().__contains__(Dimension):
            if other.__class__ is int or other.__class__ is float:
                return Derived(other / self.value, numerator=newDenom, denominator=newNum)
            else:
                raise TypeError("Can only divide float, int, other Dimension.")

        inverse = Derived(1 / self.value, numerator=newDenom, denominator=newNum)
        return other * inverse

    def __neg__(self) -> 'Derived':
        return -1*self

    def __radd__(self, other) -> 'Derived':
        return self + other

    def __rsub__(self, other) -> 'Derived':
        return -self + other

    def __copy__(self) -> 'Derived':
        return Derived(self.value,numerator=self.numerator,denominator=self.denominator)

    def __abs__(self) -> 'Derived':
        return Derived(abs(self.value),numerator=self.numerator,denominator=self.denominator)

    def __gt__(self, other) -> bool:
        #TODO: test against float
        if self.is_comparable(other):
            if not other.__class__.mro().__contains__(Dimension):
                return self.get_default() > other
            else:
                return self.get_default() > other.get_default()
        return None

    def __lt__(self, other) -> bool:
        #TODO: test against float
        if self.is_comparable(other):
            if not other.__class__.mro().__contains__(Dimension):
                return self.get_default() < other
            else:
                return self.get_default() < other.get_default()
        return None

    def __eq__(self, other) -> bool:
        #TODO: test against float
        if self.is_comparable(other):
            if not other.__class__.mro().__contains__(Dimension):
                return self.get_default() == other
            else:
                return self.get_default() == other.get_default()
        return None

    def is_dimless(self) -> bool:
        return False

    #TODO: better way to do math functions?
    def sin(self) -> float:
        if not self.is_dimless():
            raise Exception("Must be dimensionless to take sine.")
        return float(np.sin(self.value))

    def cos(self) -> float:
        if not self.is_dimless():
            raise Exception("Must be dimensionless to take cosine.")
        return float(np.cos(self.value))

    def tan(self) -> float:
        if not self.is_dimless():
            raise Exception("Must be dimensionless to take tangent.")
        return float(np.tan(self.value))

    def exp(self) -> float:
        if not self.is_dimless():
            raise Exception("Must be dimensionless to perform exponentiation.")
        return float(np.exp(self.value))

    def sqrt(self) -> 'Derived':
        return self ** (1/2)

    def ln(self) -> float:
        if not self.is_dimless():
            raise Exception("Must be dimensionless to take logarithm.")
        return float(np.log(self.value))


class Unit:
    def __init__(self, major: type, minor: str):
        self.major = major
        self.minor = minor

    def __str__(self):
        return self.major.__name__ + " " + self.minor

    def __eq__(self, other):
        return self.__class__ == other.__class__ and self.major is other.major and self.minor == other.minor

    def __hash__(self):
        return hash((self.major, self.minor))


class Derived(Dimension):

    def __init__(self, value: float, numerator=None, denominator=None):
        if denominator is None:
            denominator = []
        if numerator is None:
            numerator = []

        if numerator.__class__ is dict:
            nTemp = []
            for key in numerator.keys():
                for i in range(numerator[key]):
                    nTemp.insert(len(nTemp),key)
            numerator = nTemp
        if denominator.__class__ is dict:
            dTemp = []
            for key in denominator.keys():
                for i in range(denominator[key]):
                    dTemp.insert(len(dTemp),key)
            denominator = dTemp
        newNum = []
        newDenom = []

        unit:Unit
        for unit in numerator:
            if not unit.major.mro().__contains__(Derived) and unit.minor == unit.major.default:
                newNum.insert(len(newNum),unit)
                continue
            #weirdness here for sub derived units
            dim = unit.major(value,unit.minor)
            value = dim.get_default()
            for u in dim.numerator:
                if newDenom.__contains__(u):
                    newDenom.remove(u)
                else:
                    newNum.insert(len(newNum), Unit(u.major, u.major.default))
            for u in dim.denominator:
                if newNum.__contains__(u):
                    newNum.remove(u)
                else:
                    newDenom.insert(len(newNum), Unit(u.major, u.major.default))
        for unit in denominator:
            if not unit.major.mro().__contains__(Derived) and unit.minor == unit.major.default:
                if newNum.__contains__(unit):
                    newNum.remove(unit)
                else:
                    newDenom.insert(len(newDenom),unit)
                continue
            dim = unit.major(1/value, unit.minor)
            value = 1/dim.get_default()
            for u in dim.numerator:
                if newNum.__contains__(u):
                    newNum.remove(u)
                else:
                    newDenom.insert(len(newNum), Unit(u.major, u.major.default))
            for u in dim.denominator:
                if newDenom.__contains__(u):
                    newDenom.remove(u)
                else:
                    newNum.insert(len(newNum), Unit(u.major, u.major.default))

        self.value = value
        self.numerator = newNum
        self.denominator = newDenom

    def get_default(self) -> float:
        return self.value

    def get_special(self, units=None, derived_denom=None, warn=True) -> tuple:

        newDim = self.__copy__()

        if units is None:
            units = []

        for uni in units:
            if uni.major.mro().__contains__(Derived):
                if warn is True:
                    print("Warning. Converting to derived units are only supported in numerator.")
                derivedConverter = uni.major(1,uni.major.default)
                newDim /= derivedConverter
                newDim.numerator.append(uni)
        for uni in units:
            if uni.minor == uni.major.default:
                continue
            for unit in [u for u in newDim.numerator if u.major is uni.major]:
                dim = uni.major(1/newDim.value, uni.minor)
                newDim.value = 1/dim.get_default()
                newDim.numerator.remove(unit)
                newDim.numerator.append(uni)
            for unit in [u for u in newDim.denominator if u.major is uni.major]:
                dim = uni.major(newDim.value, uni.minor)
                newDim.value = dim.get_default()
                newDim.denominator.remove(unit)
                newDim.denominator.append(uni)

        return newDim.value, self.unit_to_str(newDim.numerator,newDim.denominator)

    def is_dimless(self) -> bool:
        return (self.numerator is None or len(self.numerator) == 0) and \
                    (self.denominator is None or len(self.denominator) == 0)


class Length(Dimension):


    default = "m"
    m = "m"
    cm = "cm"
    mm = "mm"
    inch = "inch"
    mil = "mil"
    ft = "ft"
    mile = "mile"
    nmile = "nmile"
    Au = "Au"
    pc = "pc"
    ly = "ly"
    km = "km"
    micro = "micro"

    def __init__(self, value, unit):
        self.invalid = False
        if unit == self.m:
            self.value = value
        elif unit == self.cm:
            self.value = cc.cm_m(value)
        elif unit == self.mm:
            self.value = cc.milli_(value)
        elif unit == self.inch:
            self.value = cc.inch_m(value)
        elif unit == self.mil:
            self.value = cc.mil_m(value)
        elif unit == self.ft:
            self.value = cc.ft_m(value)
        elif unit == self.mile:
            self.value = cc.mile_m(value)
        elif unit == self.nmile:
            self.value = cc.nmile_m(value)
        elif unit == self.Au:
            self.value = cc.Au_m(value)
        elif unit == self.pc:
            self.value = cc.pc_m(value)
        elif unit == self.ly:
            self.value = cc.ly_m(value)
        elif unit == self.km:
            self.value = cc.kilo_(value)
        elif unit == self.micro:
            self.value = cc.micro_(value)
        else:
            self.invalid = True
            raise ValueError("Invalid minor unit in Length.")
        self.numerator = [Unit(self.__class__, self.default)]
        self.denominator = []
        # super().__init__(self.value, numerator=[Unit(self.__class__, self.default)])

    # def __str__(self):
    #     if self.invalid:
    #         return "Invalid"
    #     else:
    #         return "{:e} m".format(self.value)

    def get_default(self):
        return self.value

    def get(self, unit: str) -> float:
        if unit == self.m:
            return self.value
        elif unit == self.cm:
            return self.value*100
        elif unit == self.mm:
            return self.value*1000
        elif unit == self.inch:
            return self.value*100/2.54
        elif unit == self.mil:
            return self.value*100*10**3/2.54
        elif unit == self.ft:
            return self.value/.3048
        elif unit == self.mile:
            return self.value/1609.344
        elif unit == self.nmile:
            return self.value/1852.0
        elif unit == self.Au:
            return self.value/1.496e11
        elif unit == self.pc:
            return self.value/3.09e16
        elif unit == self.ly:
            return self.value/9.46e15
        elif unit == self.km:
            return self.value*1e-3
        elif unit == self.micro:
            return self.value*1e6
        else:
            raise ValueError("Invalid minor unit in get Length.")


class Mass(Dimension):

    default = "kg"
    kg = "kg"
    oz = "oz"
    lbm = "lbm"
    slug = "slug"
    g = "g"
    amu = "amu"

    def __init__(self, value: float or int, unit: str):
        if value is not float or value is not int:
            if value.__class__.mro().__contains__(Dimension):
                value = float(value)
        self.invalid = False
        if unit == self.kg:
            self.value = value
        elif unit == self.oz:
            self.value = cc.oz_kg(value)
        elif unit == self.lbm:
            self.value = cc.lb_kg(value)
        elif unit == self.slug:
            self.value = cc.slug_kg(value)
        elif unit == self.g:
            self.value = value*10**(-3)
        elif unit == self.amu:
            self.value = value*cc.m_p
        else:
            self.invalid = True
            raise ValueError("Invalid minor unit in Mass.")
        self.numerator = [Unit(self.__class__, self.default)]
        self.denominator = []
        # print("init self val: ",self.value)
        # super().__init__(self.value, numerator=[Unit(self.__class__, self.default)])

    # def __str__(self):
    #     if self.invalid:
    #         return "Invalid"
    #     else:
    #         return "{:e} kg".format(self.value)

    def get_default(self):
        return self.value

    def get(self, unit: str) -> float:
        if unit == self.kg:
            return self.value
        elif unit == self.oz:
            return self.value/.02835
        elif unit == self.lbm:
            return self.value/.4536
        elif unit == self.slug:
            return self.value/14.59
        elif unit == self.g:
            return self.value*10**3
        elif unit == self.amu:
            return self.value/cc.m_p
        else:
            raise ValueError("Invalid minor unit in get Mass.")


class Time(Dimension):

    default = "s"
    s = "s"
    ns = "s"
    ms = "ms"
    minute = "minute"
    hr = "hr"
    day = "day"
    yr = "yr"
    ydhms = "ydhms"

    def __init__(self, value, unit):
        self.invalid = False
        if unit == self.s:
            self.value = value
        elif unit == self.ns:
            self.value = cc.nano_(value)
        elif unit == self.ms:
            self.value = cc.milli_(value)
        elif unit == self.minute:
            self.value = value*60
        elif unit == self.hr:
            self.value = value*60*60
        elif unit == self.day:
            self.value = value*60*60*24
        elif unit == self.yr:
            self.value = value*60*60*24*365.242
        elif unit == self.ydhms:
            print("Converting from ydhms array to seconds. May give unexpected Derived result.")
            self.value = cc.ydhms_sec(value[0],value[1],value[2],value[3],value[4],)
        else:
            self.invalid = True
            raise ValueError("Invalid minor unit in Time.")
        self.numerator = [Unit(self.__class__, self.default)]
        self.denominator = []
        # super().__init__(self.value, numerator=[Unit(self.__class__, self.default)])

    # def __str__(self):
    #     if self.invalid:
    #         return "Invalid"
    #     else:
    #         return "{:e} s".format(self.value)

    def get_default(self):
        return self.value

    def get(self, unit: str) -> float:
        if unit == self.s:
            return self.value
        elif unit == self.ns:
            return self.value * 10 ** 9
        elif unit == self.ms:
            return self.value * 10 ** 3
        elif unit == self.minute:
            return self.value/60
        elif unit == self.hr:
            return self.value/60/60
        elif unit == self.day:
            return self.value/60/60/24
        elif unit == self.yr:
            return self.value/60/60/365.242
        elif unit == self.ydhms:
            print("Converting to ydhms array. May give unexpected Derived result. "
                  "Considering 1 year = 365 days")
            total = self.value
            secs = total % 60
            total -= secs
            total /= 60
            mins = total % 60
            total -= mins
            total /= 60
            hours = total % 24
            total -= hours
            total /= 24
            days = total % 365
            total -= days
            years = total / 365
            return [years,days,hours,mins,secs]
        else:
            raise ValueError("Invalid minor unit in get Time.")


class TempEnergy(Derived):

    default = "J"
    J = "J"
    eV = "eV"
    K = "K"
    R = "R"
    C = "C"
    F = "F"
    calth = "calth"
    kcal = "kcal"
    erg = "erg"

    def __init__(self, value, unit):
        self.invalid = False
        if unit == self.J:
            self.value = value
        elif unit == self.eV:
            self.value = cc.eV_J(value)
        elif unit == self.K:
            self.value = value*cc.k_B
        elif unit == self.R:
            self.value = cc.R_K(value)*cc.k_B
        elif unit == self.C:
            print("Converting from Celsius, may give incorrect derived result. "
                  "Convert to Kelvin prior to passing in constructor.")
            self.value = cc.C_K(value)*cc.k_B
        elif unit == self.F:
            print("Converting from Fahrenheit, may give incorrect derived result. "
                  "Convert to Kelvin/Rankine prior to passing in constructor.")
            self.value = cc.F_K(value)*cc.k_B
        elif unit == self.calth:
            self.value = cc.calth_J(value)
        elif unit == self.kcal:
            self.value = cc.kcal_J(value)
        elif unit == self.erg:
            self.value = cc.ergs_J(value)
        else:
            self.invalid = True
            raise ValueError("Invalid minor unit in TempEnergy.")
        super().__init__(self.value,
                         numerator=[Unit(Mass, Mass.kg), Unit(Length,Length.m), Unit(Length,Length.m)],
                         denominator=[Unit(Time,Time.s),Unit(Time,Time.s)])

    def __str__(self):
        if self.invalid:
            return "Invalid"
        else:
            return "{:e} J".format(self.value)

    def get_default(self):
        return self.value

    def get(self, unit: str) -> float:
        if unit == self.J:
            return self.value
        elif unit == self.eV:
            return self.value / cc.e_charge
        elif unit == self.K:
            return self.value / cc.k_B
        elif unit == self.R:
            return (self.value / cc.k_B) / 1.8
        elif unit == self.C:
            print("Converting to Celsius, may give incorrect derived result.")
            return (self.value / cc.k_B) - 273.15
        elif unit == self.F:
            print("Converting to Fahrenheit, may give incorrect derived result. ")
            return ((self.value / cc.k_B) - 273.15) * (9/5) + 32
        elif unit == self.calth:
            return self.value/4.184
        elif unit == self.kcal:
            return self.value/4184
        elif unit == self.erg:
            return self.value*1e7
        else:
            raise ValueError("Invalid minor unit in get Time.")


class Power(Derived):

    default = "W"
    W = "W"

    def __init__(self, value, unit):
        self.invalid = False
        if unit == self.W:
            self.value = value
        else:
            self.invalid = True
            raise ValueError("Invalid minor unit in Power.")
        super().__init__(self.value,
                         numerator=[Unit(TempEnergy, TempEnergy.J)],
                         denominator=[Unit(Time,Time.s)])

    def __str__(self):
        if self.invalid:
            return "Invalid"
        else:
            return "{:e} W".format(self.value)

    def get_default(self):
        return self.value

    def get(self, unit: str) -> float:
        if unit == self.W:
            return self.value
        else:
            raise ValueError("Invalid minor string in get Power.")


class Pressure(Derived):

    default = "Pa"
    Pa = "Pa"
    atm = "atm"
    psi = "psi"
    torr = "torr"

    def __init__(self, value, unit):
        self.invalid = False
        if unit == self.Pa:
            self.value = value
        elif unit == self.atm:
            self.value = cc.atm_Pa(value)
        elif unit == self.psi:
            self.value = cc.psi_Pa(value)
        elif unit == self.torr:
            self.value = cc.torr_Pa(value)
        else:
            self.invalid = True
            raise ValueError("Invalid minor unit in Pressure.")
        super().__init__(self.value,
                         numerator=[Unit(Mass, Mass.kg),],
                         denominator=[Unit(Length, Length.m),Unit(Time,Time.s),Unit(Time,Time.s)])

    def __str__(self):
        if self.invalid:
            return "Invalid"
        else:
            return "{:e} Pa".format(self.value)

    def get_default(self):
        return self.value

    def get(self, unit: str) -> float:
        if unit == self.Pa:
            return self.value
        elif unit == self.psi:
            return self.value / 6894.76
        elif unit == self.atm:
            return self.value/101325.0
        elif unit == self.torr:
            return self.value/133.3223684211
        else:
            raise ValueError("Invalid minor unit in get Pressure.")


class Force(Derived):

    default = "N"
    N = "N"
    lbf = "lbf"
    dyn = "dyn"

    def __init__(self, value, unit):
        self.invalid = False
        if unit == self.N:
            self.value = value
        elif unit == self.lbf:
            self.value = cc.lbf_N(value)
        elif unit == self.dyn:
            self.value = cc.dyn_N(value)
        else:
            self.invalid = True
            raise ValueError("Invalid minor unit in Force.")
        super().__init__(self.value,
                         numerator=[Unit(Mass, Mass.kg),Unit(Length, Length.m)],
                         denominator=[Unit(Time,Time.s),Unit(Time,Time.s)])

    def __str__(self):
        if self.invalid:
            return "Invalid"
        else:
            return "{:e} N".format(self.value)

    def get_default(self):
        return self.value

    def get(self, unit: str):
        if unit == self.N:
            return self.value
        elif unit == self.lbf:
            return self.value / 4.448
        elif unit == self.dyn:
            return self.value * 1e5
        else:
            raise ValueError("Invalid minor string in get Force.")


class Temperature(Dimension):

    default = "K"
    K = "K"
    R = "R"
    C = "C"
    F = "F"

    def __init__(self, value, unit):
        self.invalid = False
        if unit == self.K:
            self.value = value
        elif unit == self.R:
            self.value = cc.R_K(value)
        elif unit == self.C:
            print("Converting from Celsius, may give incorrect derived result. "
                  "Convert to Kelvin prior to passing in constructor.")
            self.value = cc.C_K(value)
        elif unit == self.F:
            print("Converting from Fahrenheit, may give incorrect derived result. "
                  "Convert to Kelvin/Rankine prior to passing in constructor.")
            self.value = cc.F_K(value)
        else:
            self.invalid = True
            raise ValueError("Invalid minor unit in Temperature.")
        self.numerator = [Unit(self.__class__, self.default)]
        self.denominator = []
        # super().__init__(self.value, numerator=[Unit(Temperature, Temperature.K),])

    # def __str__(self):
    #     if self.invalid:
    #         return "Invalid"
    #     else:
    #         return "{:e} K".format(self.value)

    def get_default(self):
        return self.value

    def get(self, unit: str):
        if unit == self.K:
            return self.value
        elif unit == self.R:
            return self.value / 1.8
        elif unit == self.C:
            print("Converting to Celsius, may give incorrect derived result.")
            return (self.value / cc.k_B) - 273.15
        elif unit == self.F:
            print("Converting to Fahrenheit, may give incorrect derived result. ")
            return ((self.value / cc.k_B) - 273.15) * (9/5) + 32
        else:
            raise ValueError("Invalid minor string in get Temperature.")

    def to_energy(self) -> TempEnergy:
        return TempEnergy(float(Fun.k_B * self.value), TempEnergy.J)


class Angle(Derived):

    default = "rad"
    rad = "rad"
    deg = "deg"
    arcmin = "arcmin"
    arcsec = "arcsec"
    dms = "dms"

    def __init__(self, value, unit):
        self.invalid = False
        if unit == self.rad:
            self.value = value
        elif unit == self.deg:
            self.value = cc.deg_rad(value)
        elif unit == self.arcmin:
            self.value = cc.deg_rad(value/60)
        elif unit == self.arcsec:
            self.value = cc.deg_rad(value/3600)
        elif unit == self.dms:
            print("Converting from dms array to radians. May give unexpected Derived result.")
            self.value = cc.dms_rad(value[0],value[1],value[2],)
        else:
            self.invalid = True
            raise ValueError("Invalid minor unit in Angle.")
        super().__init__(self.value)

    def __str__(self):
        if self.invalid:
            return "Invalid"
        else:
            return "{:e} rad".format(self.value)

    def get_default(self):
        return self.value

    def get(self, unit: str) -> float:
        if unit == self.rad:
            return self.value
        elif unit == self.deg:
            return cc.rad_deg(self.value)
        elif unit == self.arcmin:
            return cc.rad_deg(self.value)*60
        elif unit == self.arcsec:
            return cc.rad_deg(self.value)*3600
        elif unit == self.dms:
            print("Converting to dms array. May give unexpected Derived result.")
            total = cc.rad_deg(self.value)
            degs = int(total)
            total -= degs
            amins = int(total*60)
            total -= amins/60
            asecs = total*3600
            return [degs,amins,asecs]
        else:
            raise ValueError("Invalid minor string in get Angle.")

    def is_dimless(self) -> bool:
        return True


class Current(Dimension):

    default = "A"
    A = "A"
    nA = "nA"

    def __init__(self, value, unit):
        self.invalid = False
        if unit == self.A:
            self.value = value
        elif unit == self.nA:
            self.value = cc.nano_(value)
        else:
            self.invalid = True
            raise ValueError("Invalid minor unit in Current.")
        self.numerator = [Unit(self.__class__, self.default)]
        self.denominator = []
        # super().__init__(self.value, numerator=[Unit(self.__class__, self.default)])

    # def __str__(self):
    #     if self.invalid:
    #         return "Invalid"
    #     else:
    #         return "{:e} A".format(self.value)

    def get_default(self):
        return self.value

    def get(self, unit: str) -> float:
        if unit == self.A:
            return self.value
        elif unit == self.nA:
            return self.value*10**9
        else:
            raise ValueError("Invalid minor string in get Current.")


class Charge(Derived):

    default = "C"
    C = "C"
    e = "e"

    def __init__(self, value, unit):
        self.invalid = False
        if unit == self.C:
            self.value = value
        elif unit == self.e:
            self.value = value*cc.e_charge
        else:
            self.invalid = True
            raise ValueError("Invalid minor unit in Charge.")
        super().__init__(self.value,
                         numerator=[Unit(Current, Current.A), Unit(Time,Time.s)])

    def __str__(self):
        if self.invalid:
            return "Invalid"
        else:
            return "{:e} C".format(self.value)

    def get_default(self):
        return self.value

    def get(self, unit: str) -> float:
        if unit == self.C:
            return self.value
        elif unit == self.e:
            return self.value/cc.e_charge
        else:
            raise ValueError("Invalid minor string in get Charge.")


class MagneticField(Derived):

    default = "T"
    T = "T"
    G = "G"

    def __init__(self, value, unit):
        self.invalid = False
        if unit == self.T:
            self.value = value
        elif unit == self.G:
            self.value = cc.Gauss_Tesla(value)
        else:
            self.invalid = True
            raise ValueError("Invalid minor unit in Magnetic Field.")
        super().__init__(self.value,
                         numerator=[Unit(Mass,Mass.kg)],
                         denominator=[Unit(Current, Current.A), Unit(Time,Time.s), Unit(Time,Time.s)])

    def __str__(self):
        if self.invalid:
            return "Invalid"
        else:
            return "{:e} T".format(self.value)

    def get_default(self):
        return self.value

    def get(self, unit: str) -> float:
        if unit == self.T:
            return self.value
        elif unit == self.G:
            return self.value*1e4
        else:
            raise ValueError("Invalid minor string in get MagneticField.")


class ElectricPotential(Derived):

    default = "V"
    V = "V"

    def __init__(self, value, unit):
        self.invalid = False
        if unit == self.V:
            self.value = value
        else:
            self.invalid = True
            raise ValueError("Invalid minor unit in Electric Potential.")
        super().__init__(self.value,
                         numerator={Unit(Mass,Mass.kg):1,Unit(Length,Length.m):2},
                         denominator={Unit(Current, Current.A):1, Unit(Time,Time.s):3})

    def __str__(self):
        if self.invalid:
            return "Invalid"
        else:
            return "{:e} V".format(self.value)

    def get_default(self):
        return self.value

    def get(self, unit: str) -> float:
        if unit == self.V:
            return self.value
        else:
            raise ValueError("Invalid minor string in get Electric Potential.")


class ElectricConductance(Derived):

    default = "S"
    S = "S"

    def __init__(self, value, unit):
        self.invalid = False
        if unit == self.S:
            self.value = value
        else:
            self.invalid = True
            raise ValueError("Invalid minor unit in Electric Potential.")
        super().__init__(self.value,
                         numerator={Unit(Time,Time.s):3,Unit(Current, Current.A):2},
                         denominator={Unit(Mass,Mass.kg):1,Unit(Length,Length.m):2})

    def __str__(self):
        if self.invalid:
            return "Invalid"
        else:
            return "{:e} S".format(self.value)

    def get_default(self):
        return self.value

    def get(self, unit: str) -> float:
        if unit == self.S:
            return self.value
        else:
            raise ValueError("Invalid minor string in get Electric Conductance.")


class Amount(Dimension):

    default = "mol"
    mol = "mol"
    kmol = "kmol"

    def __init__(self, value, unit):
        self.invalid = False
        if unit == self.mol:
            self.value = value
        elif unit == self.kmol:
            self.value = cc.kilo_(value)
        else:
            self.invalid = True
            raise ValueError("Invalid minor unit in Amount.")
        self.numerator = [Unit(self.__class__, self.default)]
        self.denominator = []
        # super().__init__(self.value, numerator=[Unit(Temperature, Temperature.K),])

    # def __str__(self):
    #     if self.invalid:
    #         return "Invalid"
    #     else:
    #         return "{:e} K".format(self.value)

    def get_default(self):
        return self.value

    def get(self, unit: str):
        if unit == self.mol:
            return self.value
        elif unit == self.kmol:
            return self.value / 1000
        else:
            raise ValueError("Invalid minor string in get Amount.")


class Volume(Derived):

    default = "m3"
    m3 = "m3"
    L = "L"
    mL = "mL"
    gal = "gal"

    def __init__(self, value, unit):
        self.invalid = False
        if unit == self.m3:
            self.value = value
        elif unit == self.L:
            self.value = cc.liter_m3(value)
        elif unit == self.mL:
            self.value = cc.liter_m3(cc.milli_(value))
        elif unit == self.gal:
            self.value = cc.gal_m3(value)
        else:
            self.invalid = True
            raise ValueError("Invalid minor unit in Volume.")
        super().__init__(self.value,
                         numerator={Unit(Length, Length.m):3})

    def __str__(self):
        if self.invalid:
            return "Invalid"
        else:
            return "{:e} m^3".format(self.value)

    def get_default(self):
        return self.value

    def get(self, unit: str) -> float:
        if unit == self.m3:
            return self.value
        elif unit == self.L:
            return self.value * 1e3
        elif unit == self.mL:
            return self.value * 1e6
        elif unit == self.gal:
            return self.value * 219.969
        else:
            raise ValueError("Invalid minor unit in get Volume.")


class Fun:

    # permitivity of free space - 8.85418782e-12 C^2/N*m^2
    epsilon_0 = Derived(8.85418782e-12,
                        numerator=[Unit(Charge,Charge.C),Unit(Charge,Charge.C),],
                        denominator=[Unit(Force,Force.N),Unit(Length,Length.m),Unit(Length,Length.m)])

    # Boltzmann's constant - 1.3807e-23 J/K
    k_B = Derived(1.3807e-23,
                  numerator=[Unit(TempEnergy,TempEnergy.J)],
                  denominator=[Unit(Temperature,Temperature.K)])

    #Charge of an electron - 1.6022e-19 C
    e_charge = Charge(1.6022e-19,Charge.C)

    # electron mass - 9.1093837015e-31 kg
    m_e = Mass(9.1093837015e-31,Mass.kg)

    # proton mass - 1.67262192369e-27 kg
    m_p = Mass(1.67262192369e-27,Mass.kg)

    #Earth radius at equator = 6378140 m
    R_E = Length(6378140,Length.m)

    #G times mass of earth - 3.986e14 m^3/s^2
    mu_E = Derived(3.986e14,
                   numerator=[Unit(Length,Length.m),Unit(Length,Length.m),Unit(Length,Length.m),],
                   denominator=[Unit(Time,Time.s),Unit(Time,Time.s),])

    # Newton's gravitational constant: 6.67430e-11 N*m^21/kg^2
    G = Derived(6.6743e-11,
                numerator={Unit(Force,Force.N):1,Unit(Length,Length.m):2},
                denominator={Unit(Mass,Mass.kg):2})

    # Mass of Earth: 5.9722e24 kg
    m_Earth = Mass(5.9722e24, Mass.kg)

    #Universal Gas Constant - 8.31446261815324 J/K*mol
    R = Derived(8.31446261815324,
                numerator=[Unit(TempEnergy,TempEnergy.J)],
                denominator=[Unit(Temperature,Temperature.K),Unit(Amount,Amount.mol)])

    #Plank's constant - 6.62607015 * 10^-34 J*s
    h = Derived(6.62607015e-34,
                numerator=[Unit(TempEnergy,TempEnergy.J),Unit(Time,Time.s)])

    #Speed of light in vacuum - 299792458 m/s
    c = Derived(299792458,
                numerator=[Unit(Length, Length.m)],
                denominator=[Unit(Time, Time.s)])

    #Stefan-Boltzmann constant - 2*(np.pi**5)*(k_B**4)/(15*(h**3)*(c**2))
    #TODO: something slightly wrong in reduction algorithm, had time in numerator and denominator if unit order in numerator swapped
    sigma_SB = Derived(5.6703744191844294539709967318892308758401229702913e-8,
                       numerator=[Unit(Time, Time.s), Unit(TempEnergy, TempEnergy.J), ],
                       denominator={Unit(Length, Length.m):2,Unit(Temperature, Temperature.K):4})
    #Acceleration due to gravity at Earth's surface - 9.80665 m/s^2
    g_E_0 = Derived(9.80665,
                    numerator=[Unit(Length, Length.m)],
                    denominator={Unit(Time, Time.s):2})
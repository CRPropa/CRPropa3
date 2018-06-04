from numpy import pi as M_PI

#SI units
meter = 1;
second = 1;
kilogram = 1;
ampere = 1;
mol = 1;
kelvin = 1;

#derived units
newton = 1 * kilogram * meter / second / second;
pascal = 1 * newton / meter / meter;
joule = 1 * newton * meter;
tesla = 1 * newton / ampere / meter;
volt = 1 * kilogram * meter * meter / ampere / second / second / second;
coulomb = 1 * ampere * second;

#physical constants
eplus = 1.602176487e-19 * ampere * second;
c_light = 2.99792458e8 * meter / second;
c_squared = c_light * c_light;
amu = 1.660538921e-27 * kilogram;
mass_proton = 1.67262158e-27 * kilogram;
mass_neutron = 1.67492735e-27 * kilogram;
mass_electron = 9.10938291e-31 * kilogram;
h_planck = 6.62606957e-34 * joule * second;
k_boltzmann = 1.3806488e-23 * joule / kelvin;
mu0 = 4 * M_PI * 1e-7 * newton / ampere / ampere;
epsilon0 = 1.0 / mu0 / c_squared * ampere * second / volt / meter;

#gauss
gauss = 1e-4 * tesla;
microgauss = 1e-6 * gauss;
nanogauss = 1e-9 * gauss;
muG = microgauss;
nG = nanogauss;

#electron volt
electronvolt = eplus * joule;
kiloelectronvolt = 1e3 * electronvolt;
megaelectronvolt = 1e6 * electronvolt;
gigaelectronvolt = 1e9 * electronvolt;
teraelectronvolt = 1e12 * electronvolt;
petaelectronvolt = 1e15 * electronvolt;
exaelectronvolt = 1e18 * electronvolt;
eV = electronvolt;
keV = kiloelectronvolt;
MeV = megaelectronvolt;
GeV = gigaelectronvolt;
TeV = teraelectronvolt;
PeV = petaelectronvolt;
EeV = exaelectronvolt;

#astronomical distances
au = 149597870700 * meter;
ly = 365.25 * 24 * 3600 * second * c_light;
parsec = 648000 / M_PI * au;
kiloparsec = 1e3 * parsec;
megaparsec = 1e6 * parsec;
gigaparsec = 1e9 * parsec;
pc = parsec;
kpc = kiloparsec;
Mpc = megaparsec;
Gpc = gigaparsec;

#meter
kilometer = 1000 * meter;
centimeter = 0.01 * meter;
km = kilometer;
cm = centimeter;

#second
nanosecond = 1e-9 * second;
microsecond = 1e-6 * second;
millisecond = 1e-3 * second;
minute = 60 * second;
hour = 3600 * second;
ns = nanosecond;
mus = microsecond;
ms = millisecond;
sec = second;
#ifndef UNITS_H_
#define UNITS_H_

namespace mpc {

// SI units
static const double meter = 1;
static const double second = 1;
static const double kilogram = 1;
static const double ampere = 1;
static const double mol = 1;
static const double kelvin = 1;

// derived units
static const double newton = 1;
static const double pascal = 1;
static const double joule = 1;
static const double tesla = 1;
static const double volt = 1;
static const double coulomb = 1;

// physical constants
static const double eplus = 1.602176487e-19 * ampere * second;
static const double c_light = 2.99792458e+8 * meter / second;
static const double c_squared = c_light * c_light;
static const double amu = 1.660538921e-27 * kilogram;
static const double mass_proton = 1.67262158e-27 * kilogram;

// other units

// gauss
static const double gauss = 1.e-4 * tesla;
static const double nanogauss = 1.e-9 * gauss;
static const double nG = nanogauss;

// electron volt
static const double electronvolt = eplus * joule;
static const double megaelectronvolt = 1.e6 * electronvolt;
static const double kiloelectronvolt = 1.e+3 * electronvolt;
static const double gigaelectronvolt = 1.e+9 * electronvolt;
static const double teraelectronvolt = 1.e+12 * electronvolt;
static const double petaelectronvolt = 1.e+15 * electronvolt;
static const double exaelectronvolt = 1.e+18 * electronvolt;
static const double eV = electronvolt;
static const double keV = kiloelectronvolt;
static const double MeV = megaelectronvolt;
static const double GeV = gigaelectronvolt;
static const double TeV = teraelectronvolt;
static const double PeV = petaelectronvolt;
static const double EeV = exaelectronvolt;

// parsec
static const double parsec = 3.0856775807e+16 * meter;
static const double kiloparsec = 1.e3 * parsec;
static const double megaparsec = 1.e6 * parsec;
static const double kpc = kiloparsec;
static const double Mpc = megaparsec;

} // namespace mpc

#endif /* UNITS_H_ */

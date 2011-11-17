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
static const double joule = 1;

// physical constants
static const double eplus = 1.602176487e-19 * ampere * second;
static const double c_light = 2.99792458e+8 * meter / second;
static const double c_squared = c_light * c_light;
static const double amu = 1.660538921e-27 * kilogram;

// other units
static const double electronvolt = eplus * joule;
static const double EeV = 1.e18 * electronvolt;
static const double parsec = 3.0856775807e+16 * meter;
static const double Mpc = 1.e6 * parsec;

} // namespace mpc

#endif /* UNITS_H_ */

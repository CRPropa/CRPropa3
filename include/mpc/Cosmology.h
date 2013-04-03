#ifndef MPC_COSMOLOGY_H_
#define MPC_COSMOLOGY_H_

namespace mpc {

/**
 Set the cosmological parameters for a flat universe.
 */
void setCosmologyParameters(double hubbleParameter, double omegaMatter,
		double omegaVacuum);

/**
 Hubble rate at given redshift
 H(z) = H0 * sqrt(omegaM * (1 + z)^3 + omegaL)
 */
double hubbleRate(double redshift = 0);

/**
 Redshift of a comoving object at a given comoving distance to an observer at z = 0.
 d_comoving(z) = c/H0 * int_0^z dz' / E(z')
 */
double comovingDistance2Redshift(double distance);

/**
 Comoving distance between an observer at z = 0 and a comoving object at z.
 d_comoving(z) = c/H0 * int_0^z dz' / E(z')
 */
double redshift2ComovingDistance(double redshift);

/**
 Redshift of a comoving object at a given luminosity distance to an observer at z = 0.
 d_luminosity(z) = (1 + z) * d_comoving(z)
 */
double luminosityDistance2Redshift(double distance);

/**
 Luminosity distance between an observer at z = 0 and a comoving object at z.
 d_luminosity(z) = (1 + z) * d_comoving(z)
 */
double redshift2LuminosityDistance(double redshift);

/**
 Redshift of a comoving object at a given light travel distance to an observer at z = 0.
 d_lighttravel(z) = c/H0 * int_0^z dz' / ((1 + z')  *  E(z'))
 */
double lightTravelDistance2Redshift(double distance);

/**
 Light travel distance between an observer at z = 0 and a comoving object at z.
 d_lighttravel(z) = c/H0 * int_0^z dz' / ((1 + z')  *  E(z'))
 */
double redshift2LightTravelDistance(double redshift);

} // namespace mpc

#endif // MPC_COSMOLOGY_H_

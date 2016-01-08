#include "crpropa/Cosmology.h"
#include "crpropa/Units.h"
#include "crpropa/Common.h"

#include <vector>
#include <math.h>
#include <stdexcept>

namespace crpropa {

/**
 @class Cosmology
 @brief Cosmology calculations
 */
struct Cosmology {
	double H0; // Hubble parameter at z=0
	double omegaM; // matter density parameter
	double omegaL; // vacuum energy parameter

	static const int n;
	static const double zmin;
	static const double zmax;

	std::vector<double> Z;  // redshift
	std::vector<double> Dc; // comoving distance [m]
	std::vector<double> Dl; // luminosity distance [m]
	std::vector<double> Dt; // light travel distance [m]

	void update() {
		double dH = c_light / H0; // Hubble distance

		std::vector<double> E(n);
		E[0] = 1;

		// Relation between comoving distance r and redshift z (cf. J.A. Peacock, Cosmological physics, p. 89 eq. 3.76)
		// dr = c / H(z) dz, integration using midpoint rule
		double dlz = log10(zmax) - log10(zmin);
		double dz;
		for (int i = 1; i < n; i++) {
			Z[i] = zmin * pow(10, i * dlz / (n - 1)); // logarithmic even spacing
			dz = (Z[i] - Z[i - 1]); // redshift step
			E[i] = sqrt(omegaL + omegaM * pow(1 + Z[i], 3));
			Dc[i] = Dc[i - 1] + dH * dz * (1 / E[i] + 1 / E[i - 1]) / 2;
			Dl[i] = (1 + Z[i]) * Dc[i];
			Dt[i] = Dt[i - 1]
					+ dH * dz
							* (1 / ((1 + Z[i]) * E[i])
									+ 1 / ((1 + Z[i - 1]) * E[i - 1])) / 2;
		}
	}

	Cosmology() {
        // Cosmological parameters (K.A. Olive et al. (Particle Data Group), Chin. Phys. C, 38, 090001 (2014))
		H0 = 67.3 * 1000 * meter / second / Mpc; // default values
		omegaM = 0.315;
		omegaL = 1 - omegaM;

		Z.resize(n);
		Dc.resize(n);
		Dl.resize(n);
		Dt.resize(n);

		Z[0] = 0;
		Dc[0] = 0;
		Dl[0] = 0;
		Dt[0] = 0;

		update();
	}

	void setParameters(double h, double oM) {
		H0 = h * 1e5 / Mpc;
		omegaM = oM;
		omegaL = 1 - oM;
		update();
	}
};

const int Cosmology::n = 1000;
const double Cosmology::zmin = 0.0001;
const double Cosmology::zmax = 100;

static Cosmology cosmology; // instance is created at runtime

void setCosmologyParameters(double h, double oM) {
	cosmology.setParameters(h, oM);
}

double hubbleRate(double z) {
	return cosmology.H0
			* sqrt(cosmology.omegaL + cosmology.omegaM * pow(1 + z, 3));
}

double omegaL() {
	return cosmology.omegaL;
}

double omegaM() {
	return cosmology.omegaM;
}

double H0() {
	return cosmology.H0;
}

double comovingDistance2Redshift(double d) {
	if (d < 0)
		throw std::runtime_error("Cosmology: d < 0");
	if (d > cosmology.Dc.back())
		throw std::runtime_error("Cosmology: d > dmax");
	return interpolate(d, cosmology.Dc, cosmology.Z);
}

double redshift2ComovingDistance(double z) {
	if (z < 0)
		throw std::runtime_error("Cosmology: z < 0");
	if (z > cosmology.zmax)
		throw std::runtime_error("Cosmology: z > zmax");
	return interpolate(z, cosmology.Z, cosmology.Dc);
}

double luminosityDistance2Redshift(double d) {
	if (d < 0)
		throw std::runtime_error("Cosmology: d < 0");
	if (d > cosmology.Dl.back())
		throw std::runtime_error("Cosmology: d > dmax");
	return interpolate(d, cosmology.Dl, cosmology.Z);
}

double redshift2LuminosityDistance(double z) {
	if (z < 0)
		throw std::runtime_error("Cosmology: z < 0");
	if (z > cosmology.zmax)
		throw std::runtime_error("Cosmology: z > zmax");
	return interpolate(z, cosmology.Z, cosmology.Dl);
}

double lightTravelDistance2Redshift(double d) {
	if (d < 0)
		throw std::runtime_error("Cosmology: d < 0");
	if (d > cosmology.Dt.back())
		throw std::runtime_error("Cosmology: d > dmax");
	return interpolate(d, cosmology.Dt, cosmology.Z);
}

double redshift2LightTravelDistance(double z) {
	if (z < 0)
		throw std::runtime_error("Cosmology: z < 0");
	if (z > cosmology.zmax)
		throw std::runtime_error("Cosmology: z > zmax");
	return interpolate(z, cosmology.Z, cosmology.Dt);
}

double comoving2LightTravelDistance(double d) {
	if (d < 0)
		throw std::runtime_error("Cosmology: d < 0");
	if (d > cosmology.Dc.back())
		throw std::runtime_error("Cosmology: d > dmax");
	return interpolate(d, cosmology.Dc, cosmology.Dt);
}

double lightTravel2ComovingDistance(double d) {
	if (d < 0)
		throw std::runtime_error("Cosmology: d < 0");
	if (d > cosmology.Dt.back())
		throw std::runtime_error("Cosmology: d > dmax");
	return interpolate(d, cosmology.Dt, cosmology.Dc);
}

} // namespace crpropa

#ifndef CRPROPA_HADRONICINTERACTION_H
#define CRPROPA_HADRONICINTERACTION_H

#include "crpropa/Module.h"
#include "crpropa/Vector3.h"
#include <crpropa/Grid.h>

namespace crpropa {
/**
 * \addtogroup EnergyLosses
 * @{
 */

/**
 @class HadronicInteraction
 @brief interactions of nuclei with background nucleons (Hydrogen only).
 */
class HadronicInteraction: public Module {
protected:
	double massDensity;
	ScalarGrid4d spaceTimeGrid;
	ScalarGrid spaceGrid;
	bool haveElectrons;
	bool havePhotons;
	bool haveNeutrinos;

public:
	HadronicInteraction(
		double massDensity,
		bool electrons = false,
		bool photons = false,
		bool neutrinos = false);

	HadronicInteraction(
		double massDensity,
		ScalarGrid4d spaceTimeGrid,
		bool electrons = false,
		bool photons = false,
		bool neutrinos = false);

	HadronicInteraction(
		double massDensity,
		ScalarGrid spaceGrid,
		bool electrons = false,
		bool photons = false,
		bool neutrinos = false);

	void setMassDensity(double dens);
	void setHaveElectrons(bool b);
	void setHavePhotons(bool b);
	void setHaveNeutrinos(bool b);
	void process(Candidate *candidate) const;
	double distribution_e(double energy, double x) const;
	double distribution_my1(double energy, double x) const; 
	double distribution_gamma(double energy, double x) const; 
	int numberOfElectrons(double energy) const;
	int numberOfMuonNeutrinos(double energy) const;
	int numberOfGammaRays(double energy) const;
	double CrossSection_Kelner(double energy) const;

	// these functions are not being used in the simulation
	double distribution_Carceller(double energy, double x, double jcap, double a0, double b0) const;
	double distribution_Carceller_g(double energy, double x, double jcap, double a0, double b0) const;
	double CrossSection_Carceller(double energy) const;
	double CrossSection_Galprop(double energy) const;
    Vector3d getPosition(double height, double radius) const;
};

} // namespace crpropa

#endif // CRPROPA_HADRONICINTERACTION_H

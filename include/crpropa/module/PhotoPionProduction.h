#ifndef CRPROPA_PHOTOPIONPRODUCTION_H
#define CRPROPA_PHOTOPIONPRODUCTION_H

#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"

#include <vector>

namespace crpropa {
/**
 * \addtogroup EnergyLosses
 * @{
 */

struct SophiaEventOutput {
	int nParticles;
	std::vector<double> energy;
	std::vector<int> id;
};

/**
 @class PhotoPionProduction
 @brief Photo-pion interactions of nuclei with background photons.
 */
class PhotoPionProduction: public Module {
public:
	PhotonFieldSampling photonFieldSampling;

protected:
	ref_ptr<PhotonField> photonField;
	std::vector<double> tabLorentz; ///< Lorentz factor of nucleus
	std::vector<double> tabRedshifts;  ///< redshifts (optional for haveRedshiftDependence)
	std::vector<double> tabProtonRate; ///< interaction rate in [1/m] for protons
	std::vector<double> tabNeutronRate; ///< interaction rate in [1/m] for neutrons
	double limit; ///< fraction of mean free path to limit the next step
	bool havePhotons;
	bool haveNeutrinos;
	bool haveElectrons;
	bool haveAntiNucleons;
	bool haveRedshiftDependence;

public:
	PhotoPionProduction(
		ref_ptr<PhotonField> photonField,
		bool photons = false,
		bool neutrinos = false,
		bool electrons = false,
		bool antiNucleons = false,
		double limit = 0.1,
		bool haveRedshiftDependence = false);
	void setPhotonField(ref_ptr<PhotonField> photonField);
	void setHavePhotons(bool b);
	void setHaveNeutrinos(bool b);
	void setHaveElectrons(bool b);
	void setHaveAntiNucleons(bool b);
	void setHaveRedshiftDependence(bool b);
	void setLimit(double limit);
	void initRate(std::string filename);
	double nucleonMFP(double gamma, double z, bool onProton) const;
	double nucleiModification(int A, int X) const;
	void process(Candidate *candidate) const;
	void performInteraction(Candidate *candidate, bool onProton) const;

	/**
	 Calculates the loss length E dx/dE in [m].
	 This is not used in the simulation.
	 @param	id		PDG particle id
	 @param gamma	Lorentz factor of particle
	 @param z		redshift
	 */
	double lossLength(int id, double gamma, double z = 0);

	/**
	 Direct SOPHIA interface.
	 Output is an object SophiaEventOutput with two vectors "energy" and "id" each of length N (number of out-going particles).
	 The i-th component of each vector corresponds to the same particle.
	 This is not used in the simulation.
	 @param onProton	proton or neutron
	 @param Ein			energy of nucleon
	 @param eps			energy of target photon
	 */
	SophiaEventOutput sophiaEvent(bool onProton, double Ein, double eps) const;
};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_PHOTOPIONPRODUCTION_H

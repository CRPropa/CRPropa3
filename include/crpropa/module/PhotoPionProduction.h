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

/**
 @class PhotoPionProduction
 @brief Photo-pion interactions of nuclei with background photons.
 */
class PhotoPionProduction: public Module {
protected:
	PhotonField photonField;
	Photon_Field pf;
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
		PhotonField photonField = CMB,
		bool photons = false,
		bool neutrinos = false,
		bool electrons = false,
		bool antiNucleons = false,
		double limit = 0.1,
		bool haveRedshiftDependence = false);
	void setPhotonField(PhotonField photonField);
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
	 Output is a vector of length 2*N where N is the number of out-going particles.
	 First half of the vector are the particle IDs, second half are their energies.
	 This is not used in the simulation.
	 @param onProton    proton or neutron
	 @param Ein         energy of nucleon
	 @param eps         energy of target photon
	*/
	std::vector<double> sophiaEvent(bool onProton, double Ein, double eps) const;
};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_PHOTOPIONPRODUCTION_H

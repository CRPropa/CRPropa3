#ifndef CRPROPA_PHOTOPIONPRODUCTION_H
#define CRPROPA_PHOTOPIONPRODUCTION_H

#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"

#include <vector>

namespace crpropa {

/**
 @class PhotoPionProduction
 @brief Photo-pion interactions of nuclei with background photons.
 */
class PhotoPionProduction: public Module {
private:
	PhotonField photonField;
	std::vector<double> tabLorentz; ///< Lorentz factor of nucleus
	std::vector<double> tabRedshifts; ///< Redshifts for the interaction rates
	std::vector<double> tabProtonRate; ///< interaction rate in [1/m] for protons
	std::vector<double> tabNeutronRate; ///< interaction rate in [1/m] for neutrons
	double limit; ///< fraction of mean free path to limit the next step
	bool havePhotons;
	bool haveNeutrinos;
	bool haveAntiNucleons;
	bool doRedshiftDependent;

public:
	PhotoPionProduction(PhotonField photonField = CMB, bool photons = false,
			bool neutrinos = false, bool antiNucleons = false, double limit =
					0.1);
	void setPhotonField(PhotonField photonField);
	void setHavePhotons(bool b);
	void setHaveNeutrinos(bool b);
	void setHaveAntiNucleons(bool b);
        void setDoRedshiftDependent(bool b);
	void setLimit(double limit);
	void init();
	void init(std::string filename);
	void process(Candidate *candidate) const;
	void performInteraction(Candidate *candidate, int channel) const;
	double nucleiModification(int A, int X) const;

	/**
	 Calculates the loss length E dx/dE in [m].
	 This is not used in the simulation.
	 @param	id		PDG particle id
	 @param gamma	Lorentz factor of particle
	 @param z		redshift
	 */
	double lossLength(int id, double gamma, double z = 0);
};

} // namespace crpropa

#endif // CRPROPA_PHOTOPIONPRODUCTION_H

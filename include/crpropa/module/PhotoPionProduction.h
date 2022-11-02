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
	std::string interactionTag = "PPP";

	// called by: sampleEps
	// - input: s [GeV^2]
	// - output: (s-p^2) * sigma_(nucleon/gamma) [GeV^2 * mubarn]
	double functs(double s, bool onProton) const;

	// called by: sampleEps, gaussInt
	// - input: photon energy eps [eV], Ein [GeV]
	// - output: probability to encounter photon of energy eps
	double probEps(double eps, bool onProton, double Ein, double z) const;

	/** called by: sampleEps
	@param onProton	particle type: proton or neutron
	@param Ein		energy of incoming nucleon
	- output: labframe energy [eV] of least energetic photon where PPP can occur
	 */
	double epsMinInteraction(bool onProton, double Ein) const;

	/** called by: probEps, epsMinInteraction
	@param onProton	particle type: proton or neutron
	@param Ein		energy of incoming nucleon
	- output: hadron momentum [GeV/c]
	 */
	double momentum(bool onProton, double Ein) const;
	
	// called by: functs
	// - input: photon energy [eV]
	// - output: crossection of nucleon-photon-interaction [mubarn]
	double crossection(double eps, bool onProton) const;

	// called by: crossection
	// - input: photon energy [eV], threshold [eV], max [eV], unknown [no unit]
	// - output: unknown [no unit]
	double Pl(double eps, double xth, double xMax, double alpha) const;

	// called by: crossection
	// - input: photon energy [eV], threshold [eV], unknown [eV]
	// - output: unknown [no unit]
	double Ef(double eps, double epsTh, double w) const;

	// called by: crossection
	// - input: cross section [µbarn], width [GeV], mass [GeV/c^2], rest frame photon energy [GeV]
	// - output: Breit-Wigner crossection of a resonance of width Gamma
	double breitwigner(double sigma0, double gamma, double DMM, double epsPrime, bool onProton) const;

	// called by: probEps, crossection, breitwigner, functs
	// - input: is proton [bool]
	// - output: mass [Gev/c^2]
	double mass(bool onProton) const;

	// - output: [GeV^2] head-on collision 
	double sMin() const;

	bool sampleLog = true;
	double correctionFactor = 1.6; // increeses the maximum of the propability function
	

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
	void setInteractionTag(std::string tag);
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

	/**
	 SOPHIA's photon sampling method. Returns energy [J] of a photon of the photon field.
	 @param onProton	particle type: proton or neutron
	 @param E		energy of incoming nucleon [J]
	 @param z		redshift of incoming nucleon
	 */
	double sampleEps(bool onProton, double E, double z) const;
	
	/** called by: sampleEps
	@param onProton	particle type: proton or neutron
	@param Ein		energy of incoming nucleon
	@param z		redshift of incoming nucleon
	@param epsMin   minimum photon energy of field
	@param epsMax   maximum photon energy of field
	- output: maximum probability of all photons in field
	 */
	double probEpsMax(bool onProton, double Ein, double z, double epsMin, double epsMax) const;
	
	// using log or lin spacing of photons in the range between epsMin and
	// epsMax for computing the maximum probability of photons in field
	void setSampleLog(bool log);

	// given the discrete steps to compute the maximum interaction probability pEpsMax 
	// of photons in field, the real pEpsMax may lie between the descrete tested photon energies.
	// A correction factor can be set to increase pEpsMax by that factor
	void setCorrectionFactor(double factor);

	/** get functions for the parameters of the class PhotoPionProduction, similar to the set functions */
	ref_ptr<PhotonField> getPhotonField() const;
	bool getHavePhotons() const;
	bool getHaveNeutrinos() const;
	bool getHaveElectrons() const;
	bool getHaveAntiNucleons() const;
	bool getHaveRedshiftDependence() const;
	double getLimit() const;
	bool getSampleLog() const;
	double getCorrectionFactor() const;
	std::string getInteractionTag() const;
};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_PHOTOPIONPRODUCTION_H

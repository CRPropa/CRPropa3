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
 @brief Photo-pion ("GZK") interactions of nucleons and nuclei with
 background photons.

 This module simulates the photo-pion production process, in which an
 ultra-high-energy nucleon interacts with a background photon (CMB or
 EBL/IRB) in its own rest frame, producing one or more pions:
 @f[ N + \gamma \rightarrow N' + n\,\pi @f]
 along with associated secondary photons, electrons/positrons, and
 neutrinos. For nuclei, the interaction proceeds on individual bound
 nucleons; the remaining nucleus is passed on with reduced energy while
 the participating nucleon (or the whole nucleus, depending on channel)
 is regenerated as a secondary.

 Physically, this is the interaction responsible for the
 Greisen-Zatsepin-Kuzmin (GZK) effect: above the pion-production threshold
 -- set by the invariant-mass condition for exciting a nucleon to the
 Delta(1232) resonance, @f$ s = (m_N + m_\pi)^2 @f$, reached around
 E ~ 5e19 eV for protons on CMB photons -- the interaction rate rises
 sharply. For a free nucleon undergoing resonant production, the energy
 loss per interaction corresponds to the nucleon-to-Delta(1232) mass
 ratio, roughly 1 - m_N/m_Delta ~ 24%; for a bound nucleus only the single
 interacting nucleon carries this loss, diluting the whole-nucleus loss
 fraction to roughly 1/A. This produces a "GZK horizon" of order 100 Mpc,
 beyond which UHECR protons arriving with energy above the threshold are
 exponentially suppressed -- the observed steepening of the cosmic-ray
 spectrum above ~5e19 eV. Photo-pion production is also the dominant
 mechanism for producing cosmogenic ("GZK") neutrinos, since a fraction of
 the produced pions are charged and decay into neutrinos.

 For nuclei, the per-nucleon interaction rate does not simply scale with
 mass number A: light nuclei (A <= 8) follow a reduced scaling of
 0.85 * X^(2/3) -- an exponent commonly associated with surface-shadowing
 effects, where only nucleons near the nucleus' surface interact
 efficiently -- while heavier nuclei scale approximately linearly as
 0.85 * X, closer to treating the nucleus as X quasi-independent nucleons.

 Interaction kinematics (secondary particle types, energies, and
 multiplicities) are generated with the SOPHIA Monte Carlo event generator
 (Mucke et al. 2000), an external hadronic interaction code widely used
 for photo-pion production in astroparticle physics. Interaction rates
 are precomputed as a function of nucleon Lorentz factor (and, optionally,
 redshift) and tabulated for performance, rather than evaluated from the
 full cross section at every step; the target background photon energy
 actually participating in each sampled interaction is itself drawn from
 the photon field's spectrum via rejection sampling.

 @see NuclearDecay, PhotoDisintegration, ElectronPairProduction
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
	/**
	 * @brief pion production on a given target photon field
	 * 
	 * @param photonField 	target photon field
	 * @param photons 		if true, secondary photons are added to the simulation
	 * @param neutrinos 	if true, secondary neutrinos are added to the simulation
	 * @param electrons 	if true, secondary electrons are added to the simulation
	 * @param antiNucleons 	if true, secondary anti nucleons are added to the simulation
	 * @param limit 		fraction of the mean free path, to which the propagation step will be limited
	 * @param haveRedshiftDependence 	use redshift dependent tabulated loss rates; if false, the redshift scaling of the photon field will be used
	 */
	PhotoPionProduction(
		ref_ptr<PhotonField> photonField,
		bool photons = false,
		bool neutrinos = false,
		bool electrons = false,
		bool antiNucleons = false,
		double limit = 0.1,
		bool haveRedshiftDependence = false);

	/** Sets the target background photon field for the interaction.
	 @param photonField Photon field providing the target photon spectrum
	        (e.g. CMB, EBL/IRB models). Determines both the interaction
	        rate (via its precomputed rate table) and the maximum photon
	        energy available for sampling, which sets the lowest nucleon
	        energy at which the interaction can occur at all.
	 */
	void setPhotonField(ref_ptr<PhotonField> photonField);

	// decide if secondary photons are added to the simulation
	void setHavePhotons(bool b);

	/** Decides whether secondary neutrinos (nu_e, nu_mu and antiparticles,
	 from charged pion decay) generated by SOPHIA are added to the
	 simulation as propagated candidates. These are the cosmogenic ("GZK")
	 neutrinos.
	 @param b If true, secondary neutrinos are injected.
	 */
	void setHaveNeutrinos(bool b);

	/** Decides whether secondary electrons and positrons generated by
	 SOPHIA (e.g. from pi0 decay chains and direct production) are added
	 to the simulation.
	 @param b If true, secondary electrons/positrons are injected.
	 */
	void setHaveElectrons(bool b);

	/** Decides whether secondary anti-nucleons generated by SOPHIA are
	 added to the simulation. Anti-nucleon production uses charge symmetry
	 with SOPHIA's proton/neutron output, since SOPHIA itself only
	 simulates interactions on protons and neutrons.
	 @param b If true, secondary anti-nucleons are injected.
	 */
	void setHaveAntiNucleons(bool b);

	/** Decides whether the 2-dimensional (Lorentz factor and redshift)
	 tabulated interaction rate is used instead of the default 1-dimensional
	 rate with a separate photon-field redshift scaling applied. Only
	 available for photon fields that provide a redshift-dependent
	 spectrum (e.g. EBL models); has no effect otherwise. Note that this
	 treatment is not currently available consistently across all other
	 interaction modules.
	 @param b If true, use the tabulated 2D redshift-dependent rate.
	 */
	void setHaveRedshiftDependence(bool b);

	/** Limit the propagation step to a fraction of the mean free path
	 * @param limit fraction of the mean free path
	 */	
	void setLimit(double limit);

	/** set a custom interaction tag to trace back this interaction
	 * @param tag string that will be added to the candidate and output
	 */
	void setInteractionTag(std::string tag);

	void initRate(std::string filename);

	/** get the mean free path (MFP) for a single nucleon. 
	 *  To get the MFP for the full nucleus the nucleonMFP has to be divided by by the nucleiModification factor
	 * @param gamma 	Lorentz factor of the nucleon
	 * @param z 		redshift
	 * @param onProton 	true for protons, false for neutrons
	 */
	double nucleonMFP(double gamma, double z, bool onProton) const;

	/** scaling factor for mean free path of the nucleus (converting the MFP of a single nucleon)
	 * 
	 * @param A		mass number of the nucleus
	 * @param X 	charge number of the nucleus
	 */
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

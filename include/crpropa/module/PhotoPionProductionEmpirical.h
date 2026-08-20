#ifndef CRPROPA_PHOTOPIONPRODUCTIONEMPIRICAL_H
#define CRPROPA_PHOTOPIONPRODUCTIONEMPIRICAL_H

#include "crpropa/module/PhotoPionProduction.h"

#include <string>
#include <vector>

namespace crpropa {
/**
 * \addtogroup EnergyLosses
 * @{
 */

/**
 @class PhotoPionProductionEmpirical
 @brief Photo-meson interactions using the empirical model of Morejon et al. (2019).

 An alternative to PhotoPionProduction, which treats a nucleus as a bag of
 independent nucleons (rate = nucleon rate x 0.85 A, one ejected nucleon per
 interaction).  This module implements the empirical photomeson model of
 Morejon, Fedynitch, Boncioli, Biehl & Winter, JCAP 11 (2019) 007
 [arXiv:1904.07999], as implemented in the AstroPhoMes package, which differs on
 three axes:

 - the nonelastic cross section uses a "universal function" fitted to nuclear
   photoabsorption data below 1.2 GeV and a shadowing exponent alpha(eps_r)
   running from 1.0 at threshold to 0.66 above ~1 TeV;
 - the final state follows an empirical fragmentation table, giving a mean mass
   loss <dA> = 5.8 for Fe-56 instead of 1, distributed over heavy residuals plus
   light fragments;
 - pion production scales with a separate exponent alpha_pi (2/3 at threshold,
   1 above ~50 GeV), so the meson yield per interaction is suppressed by ~3x
   below 20 GeV and enhanced by up to ~4x above it.

 SOPHIA is not used.  Secondary spectra come from the SOPHIA redistribution
 tables shipped with AstroPhoMes, and the pi/K/mu decay chain is done here.

 @note The redistribution tables are INCLUSIVE: they carry mean multiplicities
 and mean spectra with no event-by-event correlations.  Sampling them as an
 exclusive event generator is the principal approximation of this module.  Its
 consequences -- multiplicity sampling that must not manufacture impossible
 states near threshold, an energy budget that must be capped, and a charge
 balance that must be repaired -- are handled explicitly and are documented at
 the corresponding methods.

 @note Unlike PhotoPionProduction, this module needs no per-photon-field data:
 the interaction rate is built at construction from the photon field itself, so
 it works with any field, including ones defined in Python.
 */
class PhotoPionProductionEmpirical : public PhotoPionProduction {
public:
	/** How secondary multiplicities are drawn from the (inclusive) mean yields. */
	enum MultiplicitySampling {
		/** floor(m) secondaries plus one more with probability m - floor(m).
		 Exact mean at minimum variance.  The default, because the tables carry no
		 multiplicity correlations: Poisson would manufacture two-pion states in the
		 Delta region, where the mean is ~0.46 and the kinematics allow one. */
		FloorBernoulli = 0,
		/** Poisson around the mean, for studies that want the extra variance. */
		PoissonSampling = 1
	};

	/** Where the products of an unstable secondary are created. */
	enum DecayMode {
		/** at the interaction vertex, as the SOPHIA-based module does */
		PromptDecay = 0,
		/** offset by beta*gamma*c*tau along the parent direction.  Purely
		 geometric: CRPropa does not propagate mesons, so there is no energy loss,
		 magnetic deflection or cooling during the flight.  Negligible for
		 extragalactic propagation, relevant for compact source geometries. */
		DisplacedDecay = 1
	};

	/** Constructor.
	 @param photonField   target photon field
	 @param photons       if true, secondary photons are added to the simulation
	 @param neutrinos     if true, secondary neutrinos are added to the simulation
	 @param electrons     if true, secondary electrons are added to the simulation
	 @param antiNucleons  accepted for interface compatibility; the ported tables
	                      contain no anti-nucleons, so this flag does nothing
	 @param limit         fraction of the mean free path limiting the next step
	 */
	PhotoPionProductionEmpirical(
		ref_ptr<PhotonField> photonField,
		bool photons = false,
		bool neutrinos = false,
		bool electrons = false,
		bool antiNucleons = false,
		double limit = 0.1);

	/** Set the target photon field and rebuild the interaction rate.
	 The model tables themselves are field independent and are not reloaded. */
	void setPhotonField(ref_ptr<PhotonField> photonField);

	void process(Candidate *candidate) const;

	/** Interface kept for the base class; forwards to the empirical event
	 generator, which does not distinguish an incoming nucleon type. */
	void performInteraction(Candidate *candidate, bool onProton) const;

	/** Superposition scaling of the base class, neutralised: the rate returned
	 by getRate() is already the rate of the whole nucleus. */
	double nucleiModification(int A, int X) const { return 1.; }

	/** Nonelastic photomeson cross section of a nucleus.
	 Replaces the base class' crossection(eps, onProton), which is the SOPHIA
	 nucleon parameterisation.
	 @param eps  photon energy in the nucleus rest frame [J]
	 @param A    mass number
	 @param Z    charge number
	 @returns cross section [m^2]
	 */
	double crossection(double eps, int A, int Z) const;

	/** Empirical pion scaling G(A, Z, eps_r): the factor by which the pion yield
	 of a nucleus differs from the sum of its nucleons.  Equals 1 for A = 1.
	 @param eps  photon energy in the nucleus rest frame [J] */
	double pionScaling(double eps, int A, int Z) const;

	/** Interaction rate of a nucleus.
	 @param A, Z  mass and charge number
	 @param gamma Lorentz factor
	 @param z     redshift
	 @returns rate [1/m] */
	double getRate(int A, int Z, double gamma, double z = 0) const;

	/** Mean free path of a nucleus [m].  Replaces the base class' nucleonMFP. */
	double nucleusMFP(int A, int Z, double gamma, double z) const;

	/** Sample the rest-frame photon energy [J] of an interaction.  The sampling
	 density is the integrand of the rate, so sampler and rate agree by
	 construction.  Replaces the base class' sampleEps(onProton, E, z). */
	double sampleEps(int A, int Z, double gamma, double z) const;

	/** Run one interaction at a given rest-frame photon energy.  Deterministic
	 entry point for tests; process() calls it with a sampled eps. */
	void performInteraction(Candidate *candidate, double eps) const;

	/** Mean mass and charge loss per interaction implied by the fragmentation
	 table.  For Fe-56 the empirical model gives <dA> = 5.8, against 1 for the
	 superposition treatment of PhotoPionProduction. */
	void getMeanFragmentBudget(int A, int Z, double &dA, double &dZ) const;

	/** Load the eps_r grid and the cross-section basis vectors. */
	void initBasis(std::string filename);
	/** Load the SOPHIA x-spectra of the secondaries. */
	void initRedistribution(std::string filename);
	/** Load the exclusive fragmentation channels. */
	void initFragments(std::string filename);

	void setMultiplicitySampling(MultiplicitySampling m) { multiplicitySampling = m; }
	MultiplicitySampling getMultiplicitySampling() const { return multiplicitySampling; }
	void setDecayMode(DecayMode m) { decayMode = m; }
	DecayMode getDecayMode() const { return decayMode; }
	void setHaveKaons(bool b) { haveKaons = b; }
	bool getHaveKaons() const { return haveKaons; }
	void setConserveCharge(bool b) { conserveCharge = b; }
	bool getConserveCharge() const { return conserveCharge; }

	/** Number of events in which the charge could not be balanced.  Nonzero only
	 because the meson multiplicities are drawn independently; see the note on
	 inclusive sampling in the class description. */
	size_t getChargeViolationCount() const { return chargeViolations; }

	/** Number of eps_r grid points of the loaded model tables. */
	size_t getNumberOfEnergyBins() const { return tabEpsR.size(); }
	/** Rest-frame photon energy [J] of grid point i, for tests and validation. */
	double getEnergyBin(size_t i) const { return tabEpsR.at(i); }

protected:
	/** Secondary species carried by the redistribution tables, in table order. */
	enum Product { GAMMA = 0, PIPLUS, PIMINUS, PIZERO, ELECTRON, POSITRON,
	               KPLUS, KMINUS, NEUTRON, PROTON, NPRODUCTS };

	/** One x = E_secondary / E_nucleon spectrum, truncated to its contiguous
	 nonzero range and stored as a cumulative for inverse-CDF sampling. */
	struct XDist {
		double multiplicity;      ///< bin-summed yield per interaction
		int first;                ///< index of the first nonzero x bin
		std::vector<float> cdf;   ///< normalised cumulative over [first, first+size)
	};

	/** One exclusive fragmentation channel: a heavy residual plus light clusters
	 that close A and Z exactly.  Weights are energy independent. */
	struct FragmentChannel {
		short Zd, Nd;             ///< residual nucleus
		short nHe4, nHe3, nH3, nH2, nP, nN;
		double weight;
	};

	// --- model tables, all on the 159-point eps_r grid of basis.txt ----------
	// The grid is NOT log-equidistant and must not be re-gridded: the pion fade
	// windows below are hardcoded to its indices.
	std::vector<double> tabEpsR;         ///< rest-frame photon energy [J]
	std::vector<double> tabEpsRGeV;      ///< the same, in GeV
	std::vector<double> tabSigmaUniv;    ///< universal function, per nucleon [m^2]
	std::vector<double> tabSigmaProton;  ///< SOPHIA proton cross section [m^2]
	std::vector<double> tabSigmaNeutron; ///< SOPHIA neutron cross section [m^2]
	std::vector<double> tabAlpha;        ///< mass scaling of the nonelastic cross section
	std::vector<double> tabAlphaPi;      ///< mass scaling of the pion cross section
	std::vector<double> tabPionSpl;      ///< data-based pion cross section [m^2]
	std::vector<double> tabFadeLow;      ///< sigmoid weight s1, grid indices 0..31
	std::vector<double> tabFadeHigh;     ///< sigmoid weight s2, grid indices 55..104
	std::vector<double> tabWeight;       ///< trapezoid weights of the eps_r grid [J]

	// --- interaction rate, indexed [A * nlg + i] ----------------------------
	// sigma_nonel is a sum of three A-separable pieces on disjoint eps_r ranges,
	// so the rate decomposes exactly as R = R_univ + Z * R_p + N * R_n, with the
	// A dependence carried by the row index and Z applied at execution.  No
	// per-isotope and no per-photon-field table is needed.
	std::vector<double> tabRateUniv;
	std::vector<double> tabRateProton;
	std::vector<double> tabRateNeutron;

	// --- photon field flux integral I(epsMin) = int n(eps)/eps^2 deps -------
	std::vector<double> tabFluxIntegral; ///< I on a log-equidistant eps grid
	double lgFluxMin, lgFluxMax;         ///< log10 of the eps range [J]

	// --- secondary spectra and fragmentation --------------------------------
	std::vector<XDist> tabRedist[2][NPRODUCTS]; ///< [0 = neutron, 1 = proton parent]
	double lgXmin, dlgX;                        ///< log-equidistant x binning
	std::vector<std::vector<FragmentChannel> > tabFragments; ///< [Z * NSTRIDE + N]

	MultiplicitySampling multiplicitySampling;
	DecayMode decayMode;
	bool haveKaons;
	bool conserveCharge;
	mutable size_t chargeViolations;

	static const double lgmin;  ///< minimum log10(Lorentz factor)
	static const double lgmax;  ///< maximum log10(Lorentz factor)
	static const size_t nlg;    ///< number of Lorentz factor steps
	static const int Amax;      ///< largest mass number covered by the rate tables

	/** Locate eps on the eps_r grid, returning the lower index and the linear
	 weight of the upper point.  Returns false if eps is outside the grid. */
	bool locateEnergy(double eps, size_t &i, double &frac) const;

	/** sigma_nonel decomposed as u + Z * p + N * n at grid point i, so that the
	 rate integral can be factorised.  Mirrors AstroPhoMes' cs_nonel exactly,
	 including its pointwise regime masking. */
	void crossectionParts(size_t i, int A, double &u, double &p, double &n) const;

	/** Tabulate I(epsMin) for the current photon field. */
	void initFluxIntegral();
	/** Interpolate the flux integral [1/(m^3 J^2)]. */
	double fluxIntegral(double eps) const;
	/** Fill the three A-indexed rate tables from the basis and the flux integral. */
	void initRate();

	/** Draw an integer count with the requested mean; see MultiplicitySampling. */
	int sampleCount(double mean) const;
	/** Draw x = E_secondary / E_nucleon from a truncated spectrum. */
	double sampleX(const XDist &d) const;

	/** Decay an unstable secondary into CRPropa candidates, honouring the
	 have* flags and the decay mode.  `offset` accumulates the displacement of
	 the chain so far. */
	void decayPiZero(Candidate *c, double E, const Vector3d &pos) const;
	void decayPiCharged(Candidate *c, double E, int charge, const Vector3d &pos) const;
	void decayKaon(Candidate *c, double E, int charge, const Vector3d &pos) const;
	void decayMuon(Candidate *c, double E, int charge, int helicity, const Vector3d &pos) const;
	/** Position at which a secondary of energy E and rest mass m decays. */
	Vector3d decayPosition(const Candidate *c, const Vector3d &pos,
	                       double E, double m, double ctau) const;
};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_PHOTOPIONPRODUCTIONEMPIRICAL_H

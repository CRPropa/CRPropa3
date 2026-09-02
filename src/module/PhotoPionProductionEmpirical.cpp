#include "crpropa/module/PhotoPionProductionEmpirical.h"
#include "crpropa/Common.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/Random.h"
#include "kiss/logger.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace crpropa {

const double PhotoPionProductionEmpirical::lgmin = 6;
const double PhotoPionProductionEmpirical::lgmax = 14;
const size_t PhotoPionProductionEmpirical::nlg = 201;
const int PhotoPionProductionEmpirical::Amax = NUCLEAR_ZMAX + NUCLEAR_NMAX;

// energy boundaries of the empirical nonelastic cross section [GeV]
static const double EPS_SCALE = 0.3;  // mass scaling sets in
static const double EPS_UNIV = 1.2;   // universal function ends

// number of sampling points of the photon field flux integral
static const size_t N_FLUX = 1000;

// meson and lepton constants [J] and rest-frame decay lengths [m]
static const double M_PI_CHARGED = 0.13957039e9 * eV;
static const double M_PI_ZERO = 0.1349768e9 * eV;
static const double M_MUON = 0.1056583755e9 * eV;
static const double M_KAON = 0.493677e9 * eV;
static const double CTAU_PI_CHARGED = 7.8045;
static const double CTAU_MUON = 658.6384;
static const double CTAU_KAON = 3.7115;
// two-body mass ratios r = (m_mu / m_parent)^2
static const double R_PI = 0.5730487;
static const double R_KAON = 0.0457926;
// K+- branching ratios, renormalised over the two modes kept (so they add up to 1). 
// The remaining three-body and semileptonic modes (~14%) are absorbed into these; 
// kaons are ~10% of the pion yield, so the residual error on the neutrino budget is
// sub-percent.  setHaveKaons(false) removes the approximation entirely.
static const double BR_KAON_MUNU = 0.6356 / (0.6356 + 0.2066);

// product id of each table slot, matching redistribution.txt
static const int PRODUCT_ID[] = {0, 2, 3, 4, 20, 21, 50, 51, 100, 101};

PhotoPionProductionEmpirical::PhotoPionProductionEmpirical(ref_ptr<PhotonField> photonField,
		bool photons, bool neutrinos, bool electrons, bool antiNucleons, double limit)
		: PhotoPionProduction(photonField, photons, neutrinos, electrons, antiNucleons, limit) {
	setInteractionTag("PPPE");
	multiplicitySampling = FloorBernoulli;
	decayMode = PromptDecay;
	haveKaons = true;
	conserveCharge = true;
	chargeViolations = 0;
	if (antiNucleons)
		KISS_LOG_WARNING << "PhotoPionProductionEmpirical: haveAntiNucleons has no effect, "
			"the ported redistribution tables contain no anti-nucleons.\n";
	initBasis(getDataPath("PhotoPionProductionEmpirical/basis.txt"));
	initRedistribution(getDataPath("PhotoPionProductionEmpirical/redistribution.txt"));
	initFragments(getDataPath("PhotoPionProductionEmpirical/fragments.txt"));
	setPhotonField(photonField);
}

void PhotoPionProductionEmpirical::setPhotonField(ref_ptr<PhotonField> field) {
	this->photonField = field;
	setDescription("PhotoPionProductionEmpirical: " + field->getFieldName());
	// the model tables are field independent; only the rate has to be rebuilt
	if (!tabEpsR.empty()) {
		initFluxIntegral();
		initRate();
	}
}

// ---------------------------------------------------------------------------
// model tables
// ---------------------------------------------------------------------------

void PhotoPionProductionEmpirical::initBasis(std::string filename) {
	std::ifstream infile(filename.c_str());
	if (not infile.good())
		throw std::runtime_error("PhotoPionProductionEmpirical: could not open " + filename
			+ "\nThe empirical model tables are not part of the CRPropa data tarball. "
			  "See doc/pages/PhotoPionProductionEmpirical.md.");

	tabEpsR.clear(); tabEpsRGeV.clear(); tabSigmaUniv.clear();
	tabSigmaProton.clear(); tabSigmaNeutron.clear();
	tabAlpha.clear(); tabAlphaPi.clear(); tabPionSpl.clear();

	std::string line;
	while (std::getline(infile, line)) {
		if (line.empty() or line[0] == '#')
			continue;
		std::stringstream s(line);
		double e, su, sp, sn, a, api, psp;
		s >> e >> su >> sp >> sn >> a >> api >> psp;
		tabEpsRGeV.push_back(e);
		tabEpsR.push_back(e * 1e9 * eV);
		// The universal function is an extrapolating spline and goes negative
		// below its fit range (0.204 GeV); clipping it to keep above 0.
		tabSigmaUniv.push_back(std::max(su, 0.) * 1e-34);   // microbarn -> m^2
		tabSigmaProton.push_back(sp * 1e-34);
		tabSigmaNeutron.push_back(sn * 1e-34);
		tabAlpha.push_back(a);
		tabAlphaPi.push_back(api);
		// The pion spline is in microbarn.
		tabPionSpl.push_back(psp * 1e-34);
	}
	infile.close();

	const size_t n = tabEpsR.size();
	if (n < 3)
		throw std::runtime_error("PhotoPionProductionEmpirical: " + filename + " is empty");

	// trapezoid weights of the (non-equidistant) eps_r grid
	tabWeight.assign(n, 0.);
	for (size_t i = 0; i < n; i++) {
		double lo = (i == 0) ? tabEpsR[0] : 0.5 * (tabEpsR[i - 1] + tabEpsR[i]);
		double hi = (i == n - 1) ? tabEpsR[n - 1] : 0.5 * (tabEpsR[i] + tabEpsR[i + 1]);
		tabWeight[i] = hi - lo;
	}

	// The pion renormalisation of the empirical model is blended in and out using
	// two sigmoids with hardcoded parameters (AstroPhoMes photomeson_models.py:358-359).
	tabFadeLow.assign(n, 1.);
	tabFadeHigh.assign(n, 1.);
	const size_t nLow = 32, hiFirst = 55, hiLast = 105;
	for (size_t i = 0; i < n; i++) {
		if (i < nLow)
			tabFadeLow[i] = 1. / (1. + exp(-(-5. + 10. * i / (nLow - 1.))));
		if (i < hiFirst)
			tabFadeHigh[i] = 0.;
		else if (i < hiLast)
			tabFadeHigh[i] = 1. / (1. + exp(-(-5. + 10. * (i - hiFirst) / (hiLast - hiFirst - 1.))));
	}
}

bool PhotoPionProductionEmpirical::locateEnergy(double eps, size_t &i, double &frac) const {
	if (eps < tabEpsR.front() or eps > tabEpsR.back())
		return false;
	size_t hi = std::upper_bound(tabEpsR.begin(), tabEpsR.end(), eps) - tabEpsR.begin();
	if (hi == 0) hi = 1;
	if (hi >= tabEpsR.size()) hi = tabEpsR.size() - 1;
	i = hi - 1;
	frac = (eps - tabEpsR[i]) / (tabEpsR[hi] - tabEpsR[i]);
	return true;
}

// ---------------------------------------------------------------------------
// cross section
// ---------------------------------------------------------------------------

void PhotoPionProductionEmpirical::crossectionParts(size_t i, int A,
		double &u, double &p, double &n) const {
	u = p = n = 0.;
	if (A <= 1) {
		// a free nucleon is the SOPHIA cross section itself, with no scaling
		p = tabSigmaProton[i];
		n = tabSigmaNeutron[i];
		return;
	}
	// Based on AstroPhoMes cs_nonel():
	// the universal function takes over below 1.2 GeV and the mass scaling applies above 0.3 GeV.
	const double e = tabEpsRGeV[i];
	if (e <= EPS_UNIV) {
		u = A * tabSigmaUniv[i];
		if (e > EPS_SCALE)
			u *= pow(double(A), tabAlpha[i] - 1.);
	} else {
		const double s = pow(double(A), tabAlpha[i] - 1.);
		p = tabSigmaProton[i] * s;
		n = tabSigmaNeutron[i] * s;
	}
}

double PhotoPionProductionEmpirical::crossection(double eps, int A, int Z) const {
	size_t i; double f;
	if (not locateEnergy(eps, i, f))
		return 0.;
	const int N = A - Z;
	double u0, p0, n0, u1, p1, n1;
	crossectionParts(i, A, u0, p0, n0);
	crossectionParts(i + 1, A, u1, p1, n1);
	double s0 = u0 + Z * p0 + N * n0;
	double s1 = u1 + Z * p1 + N * n1;
	return std::max((1 - f) * s0 + f * s1, 0.);
}

double PhotoPionProductionEmpirical::pionScaling(double eps, int A, int Z) const {
	if (A <= 1)
		return 1.;
	size_t i; double f;
	if (not locateEnergy(eps, i, f))
		return 1.;

	const size_t j = (f < 0.5) ? i : i + 1;
	const int N = A - Z;

	// sigma^SPM_pi0 is the pi0 inclusive cross section for a superposition of
	// nucleons, it's the reference to which the data-based pion cross section
	// is renormalised.
	double spmPi0 = 0.;
	if (j < tabRedist[1][PIZERO].size())
		spmPi0 += Z * tabSigmaProton[j] * tabRedist[1][PIZERO][j].multiplicity;
	if (j < tabRedist[0][PIZERO].size())
		spmPi0 += N * tabSigmaNeutron[j] * tabRedist[0][PIZERO][j].multiplicity;
	if (spmPi0 <= 0.)
		return 1.;

	// G = A^(alpha_pi - 1) * g, with g blending the plain superposition into the
	// data-based renormalisation A * pion_spl / sigma^SPM_pi0 over the two fade
	// windows.  The pion-charge dependence cancels between the numerator and the
	// alpha_pi-rescaled denominator, so one factor serves pi+, pi- and pi0.
	const double renorm = A * tabPionSpl[j] / spmPi0;
	const double s1 = tabFadeLow[j], s2 = tabFadeHigh[j];
	double g = (1 - s2) * ((1 - s1) + s1 * renorm) + s2;
	g *= pow(double(A), tabAlphaPi[j] - 1.);
	return std::max(g, 0.);
}

// ---------------------------------------------------------------------------
// secondary spectra and fragmentation channels
// ---------------------------------------------------------------------------

void PhotoPionProductionEmpirical::initRedistribution(std::string filename) {
	std::ifstream infile(filename.c_str());
	if (not infile.good())
		throw std::runtime_error("PhotoPionProductionEmpirical: could not open " + filename
			+ "\nSee doc/pages/PhotoPionProductionEmpirical.md.");

	const size_t nE = tabEpsR.size();
	for (size_t s = 0; s < 2; s++)
		for (size_t p = 0; p < NPRODUCTS; p++)
			tabRedist[s][p].assign(nE, XDist());

	// x = E_secondary / E_parent is log-equidistant in the redistribution tables
	lgXmin = -8.025157232704402;
	double lgXmax = 0.025157232704403;
	size_t nX = 160;
	dlgX = (lgXmax - lgXmin) / nX;

	std::string line;
	while (std::getline(infile, line)) {
		if (line.empty() or line[0] == '#')
			continue;
		std::stringstream s(line);
		int nucleon, pid, k, first, n;
		double mult;
		s >> nucleon >> pid >> k >> first >> n >> mult;
		int slot = -1;
		for (size_t p = 0; p < NPRODUCTS; p++)
			if (PRODUCT_ID[p] == pid)
				slot = int(p);
		if (slot < 0 or nucleon < 0 or nucleon > 1 or k < 0 or size_t(k) >= nE)
			continue;
		XDist d;
		d.multiplicity = mult;
		d.first = first;
		d.cdf.resize(n);
		for (int i = 0; i < n; i++) {
			double v;
			s >> v;
			d.cdf[i] = float(v);
		}
		tabRedist[nucleon][slot][k] = d;
	}
	infile.close();
}

void PhotoPionProductionEmpirical::initFragments(std::string filename) {
	std::ifstream infile(filename.c_str());
	if (not infile.good())
		throw std::runtime_error("PhotoPionProductionEmpirical: could not open " + filename
			+ "\nSee doc/pages/PhotoPionProductionEmpirical.md.");

	tabFragments.clear();
	tabFragments.resize((NUCLEAR_ZMAX + 1) * NUCLEAR_NSTRIDE);

	std::string line;
	while (std::getline(infile, line)) {
		if (line.empty() or line[0] == '#')
			continue;
		std::stringstream s(line);
		int Z, N, Zd, Nd, nHe4, nHe3, nH3, nH2, nP, nN;
		double w;
		s >> Z >> N >> Zd >> Nd >> nHe4 >> nHe3 >> nH3 >> nH2 >> nP >> nN >> w;
		if (Z < 0 or Z > NUCLEAR_ZMAX or N < 0 or N > NUCLEAR_NMAX)
			continue;

		// verify A and Z conservation
		const int dA = (Z + N) - (Zd + Nd);
		const int dZ = Z - Zd;
		if (4 * nHe4 + 3 * nHe3 + 3 * nH3 + 2 * nH2 + nP + nN != dA
				or 2 * nHe4 + 2 * nHe3 + nH3 + nH2 + nP != dZ)
			throw std::runtime_error("PhotoPionProductionEmpirical: channel does not "
				"conserve A and Z in " + filename + ": " + line);

		FragmentChannel c;
		c.Zd = short(Zd); c.Nd = short(Nd);
		c.nHe4 = short(nHe4); c.nHe3 = short(nHe3); c.nH3 = short(nH3);
		c.nH2 = short(nH2); c.nP = short(nP); c.nN = short(nN);
		c.weight = w;
		tabFragments[Z * NUCLEAR_NSTRIDE + N].push_back(c);
	}
	infile.close();
}

void PhotoPionProductionEmpirical::getMeanFragmentBudget(int A, int Z,
		double &dA, double &dZ) const {
	dA = dZ = 0.;
	const int N = A - Z;
	if (Z < 0 or Z > NUCLEAR_ZMAX or N < 0 or N > NUCLEAR_NMAX)
		return;
	const std::vector<FragmentChannel> &ch = tabFragments[Z * NUCLEAR_NSTRIDE + N];
	double norm = 0.;
	for (size_t i = 0; i < ch.size(); i++) {
		dA += (A - (ch[i].Zd + ch[i].Nd)) * ch[i].weight;
		dZ += (Z - ch[i].Zd) * ch[i].weight;
		norm += ch[i].weight;
	}
	if (norm > 0.) { dA /= norm; dZ /= norm; }
}

// ---------------------------------------------------------------------------
// sampling helpers
// ---------------------------------------------------------------------------

int PhotoPionProductionEmpirical::sampleCount(double mean) const {
	if (mean <= 0.)
		return 0;
	Random &random = Random::instance();
	if (multiplicitySampling == FloorBernoulli) {
		int n = int(floor(mean));
		if (random.rand() < mean - n)
			n++;
		return n;
	}
	// Poisson sampling: 
	// using Knuth's method for small means, where its cost is acceptable, and a normal
	// approximation for large means. 
	if (mean < 30.) {
		const double L = exp(-mean);
		int k = 0;
		double p = 1.;
		do { k++; p *= random.rand(); } while (p > L);
		return k - 1;
	}
	// Note Random::randNorm's second parameter is the standard deviation despite being 
	// named "variance", so Poisson's sigma has to be passed as sqrt(mean), not mean.
	return std::max(0, int(lround(random.randNorm(mean, std::sqrt(mean)))));
}

double PhotoPionProductionEmpirical::sampleX(const XDist &d) const {
	if (d.cdf.empty())
		return 0.;
	Random &random = Random::instance();
	const size_t j = random.randBin(d.cdf);
	// log-uniform inside the bin
	return pow(10, lgXmin + (d.first + j + random.rand()) * dlgX);
}

double PhotoPionProductionEmpirical::sampleEps(int A, int Z, double gamma, double z) const {
	const size_t nE = tabEpsR.size();
	const double g = gamma * (1 + z);
	std::vector<double> cdf(nE, 0.);
	double sum = 0.;
	for (size_t i = 0; i < nE; i++) {
		double u, p, n;
		crossectionParts(i, A, u, p, n);
		const double sigma = std::max(u + Z * p + (A - Z) * n, 0.);

		sum += tabWeight[i] * tabEpsR[i] * sigma * fluxIntegral(tabEpsR[i] / (2. * g));
		cdf[i] = sum;
	}
	if (sum <= 0.)
		return 0.;
	Random &random = Random::instance();
	const size_t i = random.randBin(cdf);

	const double lo = (i == 0) ? tabEpsR[0] : 0.5 * (tabEpsR[i - 1] + tabEpsR[i]);
	const double hi = (i == nE - 1) ? tabEpsR[nE - 1] : 0.5 * (tabEpsR[i] + tabEpsR[i + 1]);
	return lo * pow(hi / lo, random.rand());
}

// ---------------------------------------------------------------------------
// photon field flux integral and interaction rate
// ---------------------------------------------------------------------------

void PhotoPionProductionEmpirical::initFluxIntegral() {
	const double eMin = photonField->getMinimumPhotonEnergy(0);
	const double eMax = photonField->getMaximumPhotonEnergy(0);
	lgFluxMin = log10(eMin);
	lgFluxMax = log10(eMax);

	// getPhotonDensity returns dn/dln(eps) [1/m^3], so dn/deps = density/eps and
	// the integrand of I(epsMin) = int dn/deps / eps^2 deps is density / eps^3.
	ref_ptr<PhotonField> field = photonField;
	auto integrand = [&field](double e) {
		return field->getPhotonDensity(e, 0.) / (e * e * e);
	};

	tabFluxIntegral.assign(N_FLUX, 0.);
	const double dlg = (lgFluxMax - lgFluxMin) / (N_FLUX - 1);
	double acc = 0.;
	for (size_t i = N_FLUX - 1; i-- > 0;) {
		double a = pow(10, lgFluxMin + dlg * i);
		double b = pow(10, lgFluxMin + dlg * (i + 1));
		acc += gaussInt(integrand, a, b);
		tabFluxIntegral[i] = acc;
	}
}

double PhotoPionProductionEmpirical::fluxIntegral(double eps) const {
	if (eps <= 0.)
		return tabFluxIntegral.front();
	const double lg = log10(eps);
	if (lg <= lgFluxMin)
		return tabFluxIntegral.front();
	if (lg >= lgFluxMax)
		return 0.;
	const double dlg = (lgFluxMax - lgFluxMin) / (N_FLUX - 1);
	const double x = (lg - lgFluxMin) / dlg;
	const size_t i = std::min(size_t(x), N_FLUX - 2);
	const double f = x - i;
	const double a = tabFluxIntegral[i], b = tabFluxIntegral[i + 1];
	// flux integral falls steeply, hence logarithmic interpolation
	if (a > 0. and b > 0.)
		return exp((1 - f) * log(a) + f * log(b));
	return (1 - f) * a + f * b;
}

void PhotoPionProductionEmpirical::initRate() {
	const size_t nE = tabEpsR.size();
	tabRateUniv.assign((Amax + 1) * nlg, 0.);
	tabRateProton.assign((Amax + 1) * nlg, 0.);
	tabRateNeutron.assign((Amax + 1) * nlg, 0.);

	// I(eps_r / 2 gamma) does not depend on A, so it is evaluated once per
	// (eps_r, gamma) pair and reused across all mass numbers.
	std::vector<double> kernel(nlg * nE);
	const double dlg = (lgmax - lgmin) / (nlg - 1);
	for (size_t j = 0; j < nlg; j++) {
		const double gamma = pow(10, lgmin + dlg * j);
		const double pre = 1. / (2. * gamma * gamma);
		for (size_t i = 0; i < nE; i++)
			kernel[j * nE + i] = pre * tabWeight[i] * tabEpsR[i]
				* fluxIntegral(tabEpsR[i] / (2. * gamma));
	}

	for (int A = 1; A <= Amax; A++) {
		for (size_t i = 0; i < nE; i++) {
			double u, p, n;
			crossectionParts(i, A, u, p, n);
			if (u == 0. and p == 0. and n == 0.)
				continue;
			for (size_t j = 0; j < nlg; j++) {
				const double k = kernel[j * nE + i];
				if (k == 0.)
					continue;
				tabRateUniv[A * nlg + j] += u * k;
				tabRateProton[A * nlg + j] += p * k;
				tabRateNeutron[A * nlg + j] += n * k;
			}
		}
	}
}

double PhotoPionProductionEmpirical::getRate(int A, int Z, double gamma, double z) const {
	if (A < 1 or A > Amax or Z < 0 or Z > A)
		return 0.;
	const double lg = log10(gamma * (1 + z));
	if (lg <= lgmin or lg >= lgmax)
		return 0.;
	const double x = (lg - lgmin) / (lgmax - lgmin) * (nlg - 1);
	const size_t j = std::min(size_t(x), nlg - 2);
	const double f = x - j;
	const size_t o = A * nlg + j;
	const int N = A - Z;
	double r0 = tabRateUniv[o] + Z * tabRateProton[o] + N * tabRateNeutron[o];
	double r1 = tabRateUniv[o + 1] + Z * tabRateProton[o + 1] + N * tabRateNeutron[o + 1];
	double rate = (1 - f) * r0 + f * r1;
	// cosmological scaling, rate per comoving distance
	rate *= pow_integer<2>(1 + z) * photonField->getRedshiftScaling(z);
	return std::max(rate, 0.);
}

double PhotoPionProductionEmpirical::nucleusMFP(int A, int Z, double gamma, double z) const {
	const double rate = getRate(A, Z, gamma, z);
	if (rate <= 0.)
		return std::numeric_limits<double>::max();
	return 1. / rate;
}

// ---------------------------------------------------------------------------
// propagation
// ---------------------------------------------------------------------------

void PhotoPionProductionEmpirical::process(Candidate *candidate) const {
	// the loop is executed at least once, to limit the next step
	double step = candidate->getCurrentStep();
	do {
		const int id = candidate->current.getId();
		if (not isNucleus(id))
			return;

		const int A = massNumber(id);
		const int Z = chargeNumber(id);
		if (A > Amax or Z > NUCLEAR_ZMAX or (A - Z) > NUCLEAR_NMAX)
			return;

		const double z = candidate->getRedshift();
		const double rate = getRate(A, Z, candidate->current.getLorentzFactor(), z);
		if (rate <= 0.)
			return;

		Random &random = Random::instance();
		const double randDistance = -log(random.rand()) / rate;
		if (step < randDistance) {
			candidate->limitNextStep(limit / rate);
			return;
		}

		performInteraction(candidate, true);
		step -= randDistance;
	} while (step > 0);
}

void PhotoPionProductionEmpirical::performInteraction(Candidate *candidate, bool) const {
	const int A = massNumber(candidate->current.getId());
	const int Z = chargeNumber(candidate->current.getId());
	const double gamma = candidate->current.getLorentzFactor();
	const double eps = sampleEps(A, Z, gamma, candidate->getRedshift());
	if (eps > 0.)
		performInteraction(candidate, eps);
}

void PhotoPionProductionEmpirical::performInteraction(Candidate *candidate, double eps) const {
	const int id = candidate->current.getId();
	const int A = massNumber(id);
	const int Z = chargeNumber(id);
	const int N = A - Z;
	const int sign = (id > 0) ? 1 : -1;
	if (Z < 0 or Z > NUCLEAR_ZMAX or N < 0 or N > NUCLEAR_NMAX)
		return;

	size_t i; double frac;
	if (not locateEnergy(eps, i, frac))
		return;

	const size_t k = (frac < 0.5) ? i : i + 1;

	const double E = candidate->current.getEnergy();
	const double EpA = E / A;
	Random &random = Random::instance();
	const Vector3d pos = random.randomInterpolatedPosition(
		candidate->previous.getPosition(), candidate->current.getPosition());

	// --- exclusive fragmentation channel --------------------------------
	// Below handles the case nucleon emission as it doesn't appear in fragments.txt. 
	// Its charge follows is determined by the meson charge and its energy from 
	// conservation. For nucleons this code should be equivalent to what SOPHIA 
	// produces.
	const bool freeNucleon = (A == 1);
	FragmentChannel fc;
	fc.Zd = short(Z); fc.Nd = short(N);
	fc.nHe4 = fc.nHe3 = fc.nH3 = fc.nH2 = fc.nP = fc.nN = 0;
	fc.weight = 1.;
	if (not freeNucleon) {
		const std::vector<FragmentChannel> &channels = tabFragments[Z * NUCLEAR_NSTRIDE + N];
		if (channels.empty())
			return;
		double cmp = random.rand(), acc = 0.;
		size_t c = channels.size() - 1;
		for (size_t j = 0; j < channels.size(); j++) {
			acc += channels[j].weight;
			if (cmp <= acc) { c = j; break; }
		}
		fc = channels[c];
	}
	const int Ares = fc.Zd + fc.Nd;
	if (Ares < 1 or (not freeNucleon and A - Ares <= 0))
		return;
	int nP = fc.nP, nN = fc.nN;

	// --- how strongly the nucleus produces mesons ------------------------
	const double sigma = crossection(eps, A, Z);
	const double spm = Z * tabSigmaProton[k] + N * tabSigmaNeutron[k];
	if (sigma <= 0. or spm <= 0.)
		return;
	// f = mean meson yield of the nucleus divided by that of its nucleons.  It
	// suppresses below ~20 GeV and enhances above, because sigma_nonel shadows as
	// A^0.66 while sigma_pi keeps scaling as A^1.
	const double f = pionScaling(eps, A, Z) * spm / sigma;
	const double wPrompt = spm / sigma;  // direct gammas and e+-, no pion scaling

	// f > 1 means more than one nucleon participates, and the energy budget has
	// to grow with it or the cap below would undo the enhancement.  The struck
	// nucleons are ejecta of the sampled channel, so mass stays exact: only their
	// energy differs from the boost-preserving E/A.
	const int nFree = nP + nN;
	int nStruck = freeNucleon ? 0
		: std::min(std::max(1, sampleCount(std::max(1., f))), nFree);

	// which parent nucleon the x spectra are read from
	const double pProton = Z * tabSigmaProton[k] / spm;

	// --- draw the broadly distributed nucleons ----------------------------
	// The channel fixes how many protons and neutrons come out, so their charge
	// is conserved by construction; only their energy is drawn here.
	std::vector<double> wideEnergy;
	std::vector<bool> wideIsProton;
	double eWide = 0.;
	int freeLeftP = nP, freeLeftN = nN;
	for (int j = 0; j < nStruck; j++) {
		const int left = freeLeftP + freeLeftN;
		if (left <= 0)
			break;
		const bool isProton = (random.rand() * left < freeLeftP);
		if (isProton) freeLeftP--; else freeLeftN--;
		const size_t parent = (random.rand() < pProton) ? 1 : 0;
		const double x = sampleX(tabRedist[parent][isProton ? PROTON : NEUTRON][k]);
		const double e = std::min(x, 1.) * EpA;
		wideEnergy.push_back(e);
		wideIsProton.push_back(isProton);
		eWide += e;
	}
	const int nWide = int(wideEnergy.size());
	// Energy released by sampling the energy for these nucleons from the SOPHIA 
	// tables, instead of E/A.	
	// For a free nucleon the leading nucleon is the residual itself, so the whole
	// energy is available and its share is what the mesons leave behind.
	double budget = freeNucleon ? 0.99 * E : nWide * EpA - eWide;
	if (budget < 0.) budget = 0.;

	// --- mesons and prompt products --------------------------------------
	std::vector<double> ePiPlus, ePiMinus, ePiZero, eKPlus, eKMinus, eGamma, eElectron, ePositron;
	double eMeson = 0.;
	struct Draw { Product slot; double weight; std::vector<double> *out; };
	const Draw draws[] = {
		{PIPLUS, f, &ePiPlus}, {PIMINUS, f, &ePiMinus}, {PIZERO, f, &ePiZero},
		{KPLUS, haveKaons ? f : 0., &eKPlus}, {KMINUS, haveKaons ? f : 0., &eKMinus},
		{GAMMA, wPrompt, &eGamma},
	};
	for (size_t d = 0; d < sizeof(draws) / sizeof(draws[0]); d++) {
		if (draws[d].weight <= 0.)
			continue;
		// mean yield of the mother, as the nucleon mix of the parent sees it
		const double mean = draws[d].weight
			* (pProton * tabRedist[1][draws[d].slot][k].multiplicity
			   + (1 - pProton) * tabRedist[0][draws[d].slot][k].multiplicity);
		const int n = sampleCount(mean);
		for (int j = 0; j < n; j++) {
			const size_t parent = (random.rand() < pProton) ? 1 : 0;
			const double x = sampleX(tabRedist[parent][draws[d].slot][k]);
			if (x <= 0.)
				continue;
			draws[d].out->push_back(x * EpA);
			eMeson += x * EpA;
		}
	}

	// Prompt e+e- come from photon conversion and Dalitz decays, so they are
	// produced as PAIRS, but the table represents them as two inclusive spectra with
	// equal multiplicity.
	{
		const double mean = wPrompt
			* (pProton * tabRedist[1][ELECTRON][k].multiplicity
			   + (1 - pProton) * tabRedist[0][ELECTRON][k].multiplicity);
		const int nPairs = sampleCount(mean);
		for (int j = 0; j < nPairs; j++) {
			const size_t pe = (random.rand() < pProton) ? 1 : 0;
			const size_t pp = (random.rand() < pProton) ? 1 : 0;
			const double xe = sampleX(tabRedist[pe][ELECTRON][k]);
			const double xp = sampleX(tabRedist[pp][POSITRON][k]);
			if (xe <= 0. or xp <= 0.)
				continue;
			eElectron.push_back(xe * EpA);
			ePositron.push_back(xp * EpA);
			eMeson += (xe + xp) * EpA;
		}
	}

	// The spectra are inclusive and uncorrelated, so their sum can fluctuate 
	// above the budget, thus rescaling to preserve the shape and energy.
	if (eMeson > budget and eMeson > 0.) {
		const double s = budget / eMeson;
		std::vector<double> *all[] = {&ePiPlus, &ePiMinus, &ePiZero, &eKPlus,
			&eKMinus, &eGamma, &eElectron, &ePositron};
		for (size_t a = 0; a < sizeof(all) / sizeof(all[0]); a++)
			for (size_t j = 0; j < all[a]->size(); j++)
				(*all[a])[j] *= s;
		eMeson = budget;
	}

	// --- charge balance ---------------------------------------------------
	// A net meson charge has to be paid for somewhere: e.g. emitting a pi+ turns a
	// proton into a neutron. Adjusting ejected nucleons first to balance the charge,
	// then from the residual nucleus itself.
	int Zres = fc.Zd;
	int qMeson = int(ePiPlus.size()) - int(ePiMinus.size())
		+ int(eKPlus.size()) - int(eKMinus.size());
	if (conserveCharge and qMeson != 0) {
		while (qMeson > 0 and nP > 0) { nP--; nN++; qMeson--; }
		while (qMeson < 0 and nN > 0) { nN--; nP++; qMeson++; }
		// a residual outside CRPropa's tabulated range is inert: PhotoDisintegration
		// and NuclearDecay both skip nuclei beyond NUCLEAR_ZMAX/NMAX
		const int ZresMax = std::min(Ares, int(NUCLEAR_ZMAX));
		const int ZresMin = std::max(0, Ares - int(NUCLEAR_NMAX));
		while (qMeson > 0 and Zres > ZresMin) { Zres--; qMeson--; }
		while (qMeson < 0 and Zres < ZresMax) { Zres++; qMeson++; }
		// last resort: kaons contribute to qMeson but have no neutral counterpart
		// to convert into, so the nuclear system cannot always absorb it
		while (qMeson > 0 and not ePiPlus.empty()) {
			ePiZero.push_back(ePiPlus.back()); ePiPlus.pop_back(); qMeson--;
		}
		while (qMeson < 0 and not ePiMinus.empty()) {
			ePiZero.push_back(ePiMinus.back()); ePiMinus.pop_back(); qMeson++;
		}
		if (qMeson != 0)
			chargeViolations++;
	}
	// p <-> n conversions preserve the ejected nucleon count, so every wide energy
	// drawn above doesn't affect the baryon number.

	// --- emit ------------------------------------------------------------
	for (int j = 0; j < fc.nHe4; j++)
		candidate->addSecondary(sign * nucleusId(4, 2), 4 * EpA, pos, 1., interactionTag);
	for (int j = 0; j < fc.nHe3; j++)
		candidate->addSecondary(sign * nucleusId(3, 2), 3 * EpA, pos, 1., interactionTag);
	for (int j = 0; j < fc.nH3; j++)
		candidate->addSecondary(sign * nucleusId(3, 1), 3 * EpA, pos, 1., interactionTag);
	for (int j = 0; j < fc.nH2; j++)
		candidate->addSecondary(sign * nucleusId(2, 1), 2 * EpA, pos, 1., interactionTag);

	size_t wj = 0;
	double eEjecta = 0.;
	for (int j = 0; j < nP; j++) {
		const double e = (wj < wideEnergy.size()) ? wideEnergy[wj++] : EpA;
		candidate->addSecondary(sign * nucleusId(1, 1), e, pos, 1., interactionTag);
		eEjecta += e;
	}
	for (int j = 0; j < nN; j++) {
		const double e = (wj < wideEnergy.size()) ? wideEnergy[wj++] : EpA;
		candidate->addSecondary(sign * nucleusId(1, 0), e, pos, 1., interactionTag);
		eEjecta += e;
	}
	eEjecta += (fc.nHe4 * 4 + fc.nHe3 * 3 + fc.nH3 * 3 + fc.nH2 * 2) * EpA;

	for (size_t j = 0; j < ePiZero.size(); j++)
		decayPiZero(candidate, ePiZero[j], pos);
	for (size_t j = 0; j < ePiPlus.size(); j++)
		decayPiCharged(candidate, ePiPlus[j], sign, pos);
	for (size_t j = 0; j < ePiMinus.size(); j++)
		decayPiCharged(candidate, ePiMinus[j], -sign, pos);
	for (size_t j = 0; j < eKPlus.size(); j++)
		decayKaon(candidate, eKPlus[j], sign, pos);
	for (size_t j = 0; j < eKMinus.size(); j++)
		decayKaon(candidate, eKMinus[j], -sign, pos);
	if (havePhotons)
		for (size_t j = 0; j < eGamma.size(); j++)
			candidate->addSecondary(22, eGamma[j], pos, 1., interactionTag);
	if (haveElectrons) {
		for (size_t j = 0; j < eElectron.size(); j++)
			candidate->addSecondary(sign * 11, eElectron[j], pos, 1., interactionTag);
		for (size_t j = 0; j < ePositron.size(); j++)
			candidate->addSecondary(-sign * 11, ePositron[j], pos, 1., interactionTag);
	}

	// --- residual nucleus -------------------------------------------------
	candidate->created = candidate->current;
	candidate->current.setId(sign * nucleusId(Ares, Zres));
	candidate->current.setEnergy(E - eEjecta - eMeson);
}

// ---------------------------------------------------------------------------
// decay chain
//
// The ported redistribution tables carry UNDECAYED mesons, and CRPropa has no
// meson transport, so the chain is done here.  All treatments are in the
// ultrarelativistic limit (E >> m always holds: the parent nucleon carries at
// least 1e6 times its rest mass), where an isotropic rest-frame decay gives a
// flat lab-energy distribution between the kinematic limits.  Every mode closes
// energy exactly.
// ---------------------------------------------------------------------------

Vector3d PhotoPionProductionEmpirical::decayPosition(const Candidate *c,
		const Vector3d &pos, double E, double m, double ctau) const {
	if (decayMode == PromptDecay or m <= 0.)
		return pos;
	Random &random = Random::instance();
	const double d = (E / m) * ctau * random.randExponential();
	return pos + c->current.getDirection() * d;
}

void PhotoPionProductionEmpirical::decayPiZero(Candidate *c, double E,
		const Vector3d &pos) const {
	if (not havePhotons)
		return;
	// pi0 -> gamma gamma; the 1.2% Dalitz mode is folded in. c*tau is 25 nm, so
	// the decay is prompt regardless of the decay mode.
	const double u = Random::instance().rand();
	c->addSecondary(22, u * E, pos, 1., interactionTag);
	c->addSecondary(22, (1 - u) * E, pos, 1., interactionTag);
}

void PhotoPionProductionEmpirical::decayPiCharged(Candidate *c, double E,
		int charge, const Vector3d &pos) const {
	Random &random = Random::instance();
	// pi -> mu nu_mu, two body: x = E_mu / E_pi is uniform on [r, 1]
	const double x = R_PI + (1 - R_PI) * random.rand();
	const double eMuon = x * E;
	const double eNu = E - eMuon;
	const Vector3d p = decayPosition(c, pos, E, M_PI_CHARGED, CTAU_PI_CHARGED);
	if (haveNeutrinos)
		c->addSecondary(charge > 0 ? 14 : -14, eNu, p, 1., interactionTag);

	// The muon helicity is fixed by x; it is the dominant shape effect on the
	// neutrino spectra, and costs one uniform deviate.
	const double fA = R_PI * (1 - x) / ((1 - R_PI) * (1 - R_PI) * x);
	const double fB = (x - R_PI) / ((1 - R_PI) * (1 - R_PI) * x);
	int helicity = (random.rand() * (fA + fB) < fA) ? +1 : -1;
	if (charge < 0)
		helicity = -helicity;
	decayMuon(c, eMuon, charge, helicity, p);
}

void PhotoPionProductionEmpirical::decayKaon(Candidate *c, double E,
		int charge, const Vector3d &pos) const {
	Random &random = Random::instance();
	const Vector3d p = decayPosition(c, pos, E, M_KAON, CTAU_KAON);
	if (random.rand() < BR_KAON_MUNU) {
		const double x = R_KAON + (1 - R_KAON) * random.rand();
		const double eMuon = x * E;
		if (haveNeutrinos)
			c->addSecondary(charge > 0 ? 14 : -14, E - eMuon, p, 1., interactionTag);
		const double fA = R_KAON * (1 - x) / ((1 - R_KAON) * (1 - R_KAON) * x);
		const double fB = (x - R_KAON) / ((1 - R_KAON) * (1 - R_KAON) * x);
		int helicity = (random.rand() * (fA + fB) < fA) ? +1 : -1;
		if (charge < 0)
			helicity = -helicity;
		decayMuon(c, eMuon, charge, helicity, p);
	} else {
		// K -> pi pi0, two-body decay with unequal masses
		const double eStar = (M_KAON * M_KAON + M_PI_CHARGED * M_PI_CHARGED
			- M_PI_ZERO * M_PI_ZERO) / (2 * M_KAON);
		const double beta = sqrt(std::max(eStar * eStar - M_PI_CHARGED * M_PI_CHARGED, 0.)) / eStar;
		const double lo = (eStar / M_KAON) * (1 - beta);
		const double hi = (eStar / M_KAON) * (1 + beta);
		const double ePi = E * (lo + (hi - lo) * random.rand());
		decayPiCharged(c, ePi, charge, p);
		decayPiZero(c, E - ePi, p);
	}
}

void PhotoPionProductionEmpirical::decayMuon(Candidate *c, double E, int charge,
		int helicity, const Vector3d &pos) const {
	if (not haveNeutrinos and not haveElectrons)
		return;
	Random &random = Random::instance();
	const Vector3d p = decayPosition(c, pos, E, M_MUON, CTAU_MUON);

	// Michel spectra in y = E_daughter / E_mu, with the helicity term (Gaisser;
	// Lipari).  P0 and P1 are the unpolarised shapes, Q0 and Q1 the polarisation
	// terms; the electron carries the opposite sign of the muon-flavour neutrino.
	// <y> is 0.35, 0.35 and 0.30 respectively, so the three sum to one in the
	// mean but not per event.
	double y[3];
	for (int s = 0; s < 3; s++) {
		const double envelope = (s == 2) ? 4.05 : 2.05;
		while (true) {
			const double a = random.rand();
			const double a2 = a * a, a3 = a2 * a;
			double g;
			if (s == 0)        // nu_mu
				g = (5. / 3 - 3 * a2 + 4 * a3 / 3) + helicity * (-1. / 3 + 3 * a2 - 8 * a3 / 3);
			else if (s == 1)   // e
				g = (5. / 3 - 3 * a2 + 4 * a3 / 3) - helicity * (-1. / 3 + 3 * a2 - 8 * a3 / 3);
			else               // nu_e
				g = (2 - 6 * a2 + 4 * a3) + helicity * (2 - 12 * a + 18 * a2 - 8 * a3);
			if (random.rand() * envelope < g) { y[s] = a; break; }
		}
	}
	const double norm = y[0] + y[1] + y[2];
	if (norm <= 0.)
		return;
	const double eNuMu = E * y[0] / norm;
	const double eLepton = E * y[1] / norm;
	const double eNuE = E * y[2] / norm;

	if (haveNeutrinos) {
		// mu+ -> e+ nu_e anti-nu_mu;  mu- -> e- anti-nu_e nu_mu
		c->addSecondary(charge > 0 ? -14 : 14, eNuMu, p, 1., interactionTag);
		c->addSecondary(charge > 0 ? 12 : -12, eNuE, p, 1., interactionTag);
	}
	if (haveElectrons)
		c->addSecondary(charge > 0 ? -11 : 11, eLepton, p, 1., interactionTag);
}

}

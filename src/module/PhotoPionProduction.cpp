#include "crpropa/module/PhotoPionProduction.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/Random.h"
#include "crpropa/PhotonBackground.h"

#include <kiss/convert.h>
#include <kiss/logger.h>
#include "sophia.h"

#include <limits>
#include <cmath>
#include <sstream>
#include <fstream>
#include <stdexcept>

namespace crpropa {


PhotoPionProduction::PhotoPionProduction( PhotonField field,
                                          bool photons,
                                          bool neutrinos,
                                          bool electrons,
                                          bool antiNucleons,
                                          std::string tag,
                                          bool useTabData,
                                          double l) {
    setPhotonField(field);
    this->customPhotonField = CustomPhotonField(getDataPath("Scaling/" + photonFieldName(field) + ".txt"));
    this->spaceTimeGrid = ScalarGrid4d();
    this->spaceGrid = ScalarGrid();
    havePhotons = photons;
    haveNeutrinos = neutrinos;
    haveElectrons = electrons;
    haveAntiNucleons = antiNucleons;
    this-> tag = tag;
    useTabulatedData = useTabData;
    if (useTabData) initHistogram(getDataPath("PhotoPionProduction/SOPHIA_histogram.txt"));
    limit = l;
    setDescription("PhotoPionProduction_isotropicConstant");
}

PhotoPionProduction::PhotoPionProduction( PhotonField field,
                                          ScalarGrid4d spaceTimeGrid,
                                          bool photons,
                                          bool neutrinos,
                                          bool electrons,
                                          bool antiNucleons,
                                          std::string tag,
                                          bool useTabData,
                                          double l) {
    setPhotonField(field);
    this->customPhotonField = CustomPhotonField(getDataPath("Scaling/" + photonFieldName(field) + ".txt"));
    this->spaceTimeGrid = spaceTimeGrid;
    this->spaceGrid = ScalarGrid();
    havePhotons = photons;
    haveNeutrinos = neutrinos;
    haveElectrons = electrons;
    haveAntiNucleons = antiNucleons;
    this-> tag = tag;
    useTabulatedData = useTabData;
    if (useTabData) initHistogram(getDataPath("PhotoPionProduction/SOPHIA_histogram.txt"));
    limit = l;
    setDescription("PhotoPionProduction_spaceDependentConstant");
}

PhotoPionProduction::PhotoPionProduction( PhotonField field,
                                          ScalarGrid spaceGrid,
                                          bool photons,
                                          bool neutrinos,
                                          bool electrons,
                                          bool antiNucleons,
                                          std::string tag,
                                          bool useTabData,
                                          double l) {
    setPhotonField(field);
    this->customPhotonField = CustomPhotonField(getDataPath("Scaling/" + photonFieldName(field) + ".txt"));
    this->spaceTimeGrid = ScalarGrid4d();
    this->spaceGrid = spaceGrid;
    havePhotons = photons;
    haveNeutrinos = neutrinos;
    haveElectrons = electrons;
    haveAntiNucleons = antiNucleons;
    this-> tag = tag;
    useTabulatedData = useTabData;
    if (useTabData) initHistogram(getDataPath("PhotoPionProduction/SOPHIA_histogram.txt"));
    limit = l;
}

void PhotoPionProduction::setPhotonField(PhotonField field) {
    photonField = field;
    std::string fname = photonFieldName(field);
    initRate(getDataPath("PhotoPionProduction/rate_" + fname + ".txt"));
    this->customPhotonField = CustomPhotonField(getDataPath("Scaling/" + photonFieldName(field) + ".txt"));
}

void PhotoPionProduction::setHavePhotons(bool b) {
    havePhotons = b;
}

void PhotoPionProduction::setHaveElectrons(bool b) {
    haveElectrons = b;
}

void PhotoPionProduction::setHaveNeutrinos(bool b) {
    haveNeutrinos = b;
}

void PhotoPionProduction::setHaveAntiNucleons(bool b) {
    haveAntiNucleons = b;
}

void PhotoPionProduction::setUseTabulatedData(bool b) {
    useTabulatedData = b;
}

void PhotoPionProduction::setLimit(double l) {
    limit = l;
}

void PhotoPionProduction::initRate(std::string filename) {
    // clear previously loaded tables
    tabLorentz.clear();
    tabRedshifts.clear();
    tabProtonRate.clear();
    tabNeutronRate.clear();

    std::ifstream infile(filename.c_str());
    if (!infile.good())
        throw std::runtime_error("PhotoPionProduction: could not open file " + filename);

    double zOld = -1, aOld = -1;
    while (infile.good()) {
        if (infile.peek() == '#') {
            infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            continue;
        }
        double z, a, b, c;
        infile >> z >> a >> b >> c;
        if (!infile)
            break;
        if (z > zOld) {
            tabRedshifts.push_back(z);
            zOld = z;
        }
        if (a > aOld) {
            tabLorentz.push_back(pow(10, a));
            aOld = a;
        }
        tabProtonRate.push_back(b / Mpc);
        tabNeutronRate.push_back(c / Mpc);
    }

    infile.close();
}

/*
    related to histogram version of SOPHIA:
        - initHistogram
        - hashTag
        - produce
        - drawEnergy
        - snapToHalfLog
        - sophiaEvent
*/

void PhotoPionProduction::initHistogram(std::string filename) {
    // read in histogram file
    hashMap.clear();
    histData.clear();
    std::ifstream infile(filename.c_str());
    if (!infile.good())
        throw std::runtime_error("PhotoPionProduction: Could not open file " + filename);

    std::string line;
    while (std::getline(infile, line)) {
        std::istringstream ss(line);
        std::string hash;
        ss >> hash;
        std::vector<double> vec;
        double n;
        while (ss >> n)
            vec.push_back(n);
        // input.insert({hash, vec});  // with C++11
        hashMap.push_back(hash);
        histData.push_back(vec);
    }
    infile.close();
}


std::string PhotoPionProduction::hashTag(int n,     // nucleon type
                                         double E,  // primary Energy
                                         double e,  // photon energy
                                         int ID,    // particle to be produced
                                         int m = 0  // amount of particles produced in this event
                                         ) const {
    // method to generate hash tags for navigation in std::unordered_map
    std::stringstream hash;
    int exp = std::floor(log10(E));
    int pre = E / pow(10, exp);
    hash << "#" << n << "_" << pre << "e";
    (exp >= 0)? hash << "+" : hash << "-";
    if (exp < 10) hash << "0";
    hash << std::abs(exp) << "_";
    exp = std::floor(log10(e));
    pre = e / pow(10, exp);
    hash << pre << "e";
    (exp >= 0)? hash << "+" : hash << "-";
    if (exp < 10) hash << "0";
    hash << std::abs(exp) << "_" << ID;
    if (m > 0) hash << "_x" << m;
    return hash.str();
}


int PhotoPionProduction::produce(const std::vector<double> &particle) const {
/*
    - input: probability vector with chances to produce: [0,1,...n] particles
    - output: number of particles being produced
*/
    if (particle.size() == 0)
        return 0;
    Random &random = Random::instance();
    double r = random.rand();
    int index = 0;
    while ((r >= 0) && (index < particle.size())) {
        r -= particle[index];
        index++;
    }
    return index - 1;
}


double PhotoPionProduction::drawEnergy(const std::vector<double> &data) const {
    /*
        input format: first half of vector contains probabilities to draw
                      a certain energy contained in the second half
        output: energy
    */
    if (data.size() == 0)
        return 0.;
    std::vector<double> p, E;
    for (int i = 0; i < data.size(); ++i) {
        p.push_back(data[i]);
        E.push_back(data[i + data.size() / 2]);
    }
    int pos = produce(p);
    if (pos == 0)
        return E[pos];
    // interpolation
    Random &random = Random::instance();
    double r = random.rand();
    return E[pos - 1] * (1. - r) + r * E[pos];
}


double PhotoPionProduction::snapToHalfLog(double x) const {
    /*
        method to aid the hashTag function
        selects the closest value where histogram data is available
    */
    int exp = std::floor(log10(x));
    double pre = x / pow(10, exp);
    if (pre == 1.0)
        return x;
    double result = pow(10, std::ceil(log10(x)));
    if (pre < 2.5)
        return result / 10.;
    if (pre >= 7.5)
        return result;
    return result / 2.;
}


std::vector<double> PhotoPionProduction::sophiaEvent(bool onProton,  // 0=p, 1=n
                                                     double E,       // primary nucleon's energy / GeV
                                                     double e        // target photon's energy / eV
                                                     ) const {
    /*
        Histogram version of SOPHIA.
    */
    const double nature = 1 - static_cast<int>(onProton);
    const double E_in = snapToHalfLog(E);
    const double eps = snapToHalfLog(e);
    std::vector<double> output;
    std::string hash;
    ptrdiff_t whereHash;

    hash = hashTag(nature, E_in, eps, 13);
    whereHash = std::find(hashMap.begin(), hashMap.end(), hash) - hashMap.begin();
    const std::vector<double> proton = histData[whereHash];

    hash = hashTag(nature, E_in, eps, 14);
    whereHash = std::find(hashMap.begin(), hashMap.end(), hash) - hashMap.begin();
    const std::vector<double> neutron = histData[whereHash];

    // return primary if no histogram for E and e is available
    if (proton.size() == 0 && neutron.size() == 0) {
        int id = (onProton)? 13 : 14;
        output.push_back(id);
        output.push_back(E_in);
        return output;
    }

    hash = hashTag(nature, E_in, eps, 1);
    whereHash = std::find(hashMap.begin(), hashMap.end(), hash) - hashMap.begin();
    const std::vector<double> photon = histData[whereHash];
    
    hash = hashTag(nature, E_in, eps, 2);
    whereHash = std::find(hashMap.begin(), hashMap.end(), hash) - hashMap.begin();
    const std::vector<double> positron = histData[whereHash];

    hash = hashTag(nature, E_in, eps, 3);
    whereHash = std::find(hashMap.begin(), hashMap.end(), hash) - hashMap.begin();
    const std::vector<double> electron = histData[whereHash];

    hash = hashTag(nature, E_in, eps, -13);
    whereHash = std::find(hashMap.begin(), hashMap.end(), hash) - hashMap.begin();
    const std::vector<double> antiProton = histData[whereHash];

    hash = hashTag(nature, E_in, eps, -14);
    whereHash = std::find(hashMap.begin(), hashMap.end(), hash) - hashMap.begin();
    const std::vector<double> antiNeutron = histData[whereHash];

    hash = hashTag(nature, E_in, eps, 15);
    whereHash = std::find(hashMap.begin(), hashMap.end(), hash) - hashMap.begin();
    const std::vector<double> nu_e = histData[whereHash];

    hash = hashTag(nature, E_in, eps, 16);
    whereHash = std::find(hashMap.begin(), hashMap.end(), hash) - hashMap.begin();
    const std::vector<double> antiNu_e = histData[whereHash];

    hash = hashTag(nature, E_in, eps, 17);
    whereHash = std::find(hashMap.begin(), hashMap.end(), hash) - hashMap.begin();
    const std::vector<double> nu_mu = histData[whereHash];

    const std::vector<double> antiNu_mu = nu_mu;

// ########################################################
// ### particle production (without energy)
// ########################################################
    const int Le_in = 0;
    const int charge_in = 1 - nature;
    const int Nnuc_in = 1;

    int Le = 0, charge = 0, Nnuc = 0;
    int N_electron = 0, N_positron = 0,
        N_antiNu_e = 0, N_nu_e = 0,
        N_proton = 0, N_neutron = 0,
        N_antiProton = 0, N_antiNeutron = 0,
        N_photon = 0, N_nu_mu = 0, N_antiNu_mu = 0;

    // looped particle production (preserves quantum numbers)
    do {  // lepton loop
        do {  // charge loop
            do {  // nucleon loop
                // reset particles for new production loop
                Le = 0;
                charge = 0;
                Nnuc = 0;
                N_electron = 0;
                N_positron = 0;
                N_antiNu_e = 0;
                N_nu_e = 0;
                N_proton = 0;
                N_neutron = 0;
                N_antiProton = 0;
                N_antiNeutron = 0;

                N_proton += produce(proton);
                charge += N_proton;
                Nnuc += N_proton;

                N_neutron += produce(neutron);
                Nnuc += N_neutron;

                N_nu_e += produce(nu_e);
                Le += N_nu_e;

                N_antiNu_e += produce(antiNu_e);
                Le -= N_antiNu_e;

                N_electron += produce(electron);
                charge -= N_electron;
                Le += N_electron;

                N_positron += produce(positron);
                charge += N_positron;
                Le -= N_positron;

                N_antiProton += produce(antiProton);
                charge -= N_antiProton;
                Nnuc -= N_antiProton;

                N_antiNeutron += produce(antiNeutron);
                Nnuc -= N_antiNeutron;
            } while ( Nnuc != Nnuc_in );
        } while ( charge != charge_in );
    } while ( Le != Le_in );
    do {
        N_nu_mu = 0;
        N_antiNu_mu = 0;
        N_nu_mu += produce(nu_mu);
        N_antiNu_mu += N_nu_mu;  // antiNu_mu equals nu_mu
    } while ( N_nu_mu != (N_nu_e + N_antiNu_e) );

    // experimental setup
    // N_nu_mu = N_nu_e + N_antiNu_e;
    // N_antiNu_mu = N_nu_mu;  // antiNu_mu equals nu_mu
    
    N_photon = produce(photon);

    // std::cout << "[" << N_electron << "e-," << N_positron << "e+," << N_nu_e << "ne," << N_antiNu_e << "~ne," << N_proton << "p,"
    //           << N_neutron << "n," << N_nu_mu << "nm," << N_antiNu_mu << "~nm,"
    //           << N_photon << "ph," << N_antiProton << "~p," << N_antiNeutron << "~n]@["
    //           << Le << "L," << charge << "C," << Nnuc << "N] " << std::endl;

// ########################################################
// ### draw energy of produced particles
// ########################################################
    double availableEnergy = E;
    std::vector<double> outE;
    double E_part;
    const int pCount[] = {N_antiNeutron, N_antiProton, N_photon,
                          N_positron, N_electron, N_proton,
                          N_neutron, N_nu_e, N_antiNu_e,
                          N_nu_mu, N_antiNu_mu};
    const int partID[] = {-14, -13, 1, 2, 3, 13, 14, 15, 16, 17, 18};
    for (int j = 0; j < 11; ++j) {
        for (int k = 0; k < pCount[j]; ++k) {
            output.push_back(partID[j]);
            int id = (partID[j] == 18)? 17 : partID[j];  // ~nu_mu=nu_mu
            hash = hashTag(nature, E_in, eps, id, pCount[j]);
            whereHash = std::find(hashMap.begin(), hashMap.end(), hash) - hashMap.begin(); 
            E_part = drawEnergy(histData[whereHash]);
            availableEnergy -= E_part;
            outE.push_back(E_part);
        }
    }
    // preserve energy
    double weight = E / (E - availableEnergy);
    int nOutPart = output.size();
    for (int j = 0; j < nOutPart; ++j) {
        output.push_back(outE[j] * weight);
    }
    return output;
}


double PhotoPionProduction::nucleonMFP(double gamma, double z, bool onProton, Vector3d pos, double time) const {
    const std::vector<double> &tabRate = (onProton)? tabProtonRate : tabNeutronRate;

    // scale nucleus energy instead of background photon energy
    gamma *= (1 + z);
    if (gamma < tabLorentz.front() or (gamma > tabLorentz.back()))
        return std::numeric_limits<double>::max();

    // geometric scaling
    double rate = 1.;
    const std::string description = getDescription();
    if (description == "PhotoPionProduction_isotropicConstant") {
        // do nothing, just check for correct initialization
    } else if (description == "PhotoPionProduction_spaceDependentConstant") {
        rate *= spaceGrid.interpolate(pos);
    } else if (description == "PhotoPionProduction_spaceTimeDependent") {
        rate *= spaceTimeGrid.interpolate(pos, time);
    } else {
        throw std::runtime_error("PhotoPionProduction: invalid description string");
    }
    if (rate == 0.)
        return std::numeric_limits<double>::max();

    rate *= interpolate2d(z, gamma, tabRedshifts, tabLorentz, tabRate);

    // cosmological scaling
    rate *= pow(1 + z, 2);

    return 1. / rate;
}

double PhotoPionProduction::nucleiModification(int A, int X) const {
    if (A == 1)
        return 1.;
    if (A <= 8)
        return 0.85 * pow(X, 2. / 3.);
    return 0.85 * X;
}

void PhotoPionProduction::process(Candidate *candidate) const {
    double step = candidate->getCurrentStep();
    double z = candidate->getRedshift();
    Vector3d pos = candidate->current.getPosition();
    double time = candidate->getTrajectoryLength()/c_light;
    // the loop is processed at least once for limiting the next step
    do {
        // check if nucleus
        int id = candidate->current.getId();
        if (!isNucleus(id))
            return;

        // find interaction with minimum random distance
        Random &random = Random::instance();
        double randDistance = std::numeric_limits<double>::max();
        double meanFreePath;
        double totalRate = 0;
        bool onProton = true; // interacting particle: proton or neutron

        int A = massNumber(id);
        int Z = chargeNumber(id);
        int N = A - Z;
        double gamma = candidate->current.getLorentzFactor();

        // check for interaction on protons
        if (Z > 0) {
            meanFreePath = nucleonMFP(gamma, z, true, pos, time) / nucleiModification(A, Z);
            randDistance = -log(random.rand()) * meanFreePath;
            totalRate += 1. / meanFreePath;
        }
        // check for interaction on neutrons
        if (N > 0) {
            meanFreePath = nucleonMFP(gamma, z, false, pos, time) / nucleiModification(A, N);
            totalRate += 1. / meanFreePath;
            double d = -log(random.rand()) * meanFreePath;
            if (d < randDistance) {
                randDistance = d;
                onProton = false;
            }
        }
        // check if interaction does not happen
        if ( meanFreePath == std::numeric_limits<double>::max())
            return;
        if (step < randDistance) {
            if (totalRate > 0.)
                candidate->limitNextStep(limit / totalRate);
            return;
        }
        // interact and repeat with remaining step
        performInteraction(candidate, onProton);
        step -= randDistance;
    } while (step > 0);
}

void PhotoPionProduction::performInteraction(Candidate *candidate, bool onProton) const {
    int id = candidate->current.getId();
    int A = massNumber(id);
    int Z = chargeNumber(id);
    double E = candidate->current.getEnergy();
    double EpA = E / A;
    double z = candidate->getRedshift();

    // SOPHIA simulates interactions only for protons / neutrons
    // for anti-protons / neutrons assume charge symmetry and change all
    // interaction products from particle <--> anti-particle
    int sign = (id > 0) ? 1 : -1;

    // SOPHIA - input:
    int nature = 1 - static_cast<int>(onProton);  // 0=proton, 1=neutron
    double Ein = EpA / GeV;
    double eps = customPhotonField.sampleEps(onProton, Ein, z);

    // SOPHIA - output:
    double outputEnergy[2000];
    int outPartID[2000];
    int nOutPart;

    if (useTabulatedData) {
        std::vector<double> outVec = sophiaEvent(onProton, Ein, eps);
        nOutPart = outVec.size() / 2;
        for (int i = 0; i < nOutPart; ++i) {
            outPartID[i] = outVec[i];
            outputEnergy[i] = outVec[i+nOutPart];
        }
    } else {
        #pragma omp critical
        {
            sophiaevent_(nature, Ein, eps, outputEnergy, outPartID, nOutPart);
        }
    }

    // output particle treatment
	Random &random = Random::instance();
	Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
	std::vector<int> pnType;  // filled with either 13 (proton) or 14 (neutron)
	std::vector<double> pnEnergy;  // corresponding energies of proton or neutron
	for (int i = 0; i < nOutPart; i++) { // loop over out-going particles
		double Eout = outputEnergy[i] * GeV; // only the energy is used; could be changed for more detail
		int pType = outPartID[i];
		switch (pType) {
		case 13: // proton
		case 14: // neutron
			// proton and neutron data is taken to determine primary particle in a later step
			pnType.push_back(pType);
			pnEnergy.push_back(Eout);
			break;
		case -13: // anti-proton
		case -14: // anti-neutron
			if (haveAntiNucleons)
				try
				{
					candidate->addSecondary(-sign * nucleusId(1, 14 + pType), Eout, pos, tag);
				}
				catch (std::runtime_error &e)
				{
					KISS_LOG_ERROR<< "Something went wrong in the PhotoPionProduction (anti-nucleon production)\n" << "Something went wrong in the PhotoPionProduction\n"<< "Please report this error on https://github.com/CRPropa/CRPropa3/issues including your simulation setup and the following random seed:\n" << Random::instance().getSeed_base64();
					throw;
				}
			break;
		case 1: // photon
			if (havePhotons)
				candidate->addSecondary(22, Eout, pos, tag);
			break;
		case 2: // positron
			if (haveElectrons)
				candidate->addSecondary(sign * -11, Eout, pos, tag);
			break;
		case 3: // electron
			if (haveElectrons)
				candidate->addSecondary(sign * 11, Eout, pos, tag);
			break;
		case 15: // nu_e
			if (haveNeutrinos)
				candidate->addSecondary(sign * 12, Eout, pos, tag);
			break;
		case 16: // antinu_e
			if (haveNeutrinos)
				candidate->addSecondary(sign * -12, Eout, pos, tag);
			break;
		case 17: // nu_muon
			if (haveNeutrinos)
				candidate->addSecondary(sign * 14, Eout, pos, tag);
			break;
		case 18: // antinu_muon
			if (haveNeutrinos)
				candidate->addSecondary(sign * -14, Eout, pos, tag);
			break;
		default:
			throw std::runtime_error("PhotoPionProduction: unexpected particle " + kiss::str(pType));
		}
	}

    // threshold check is removed from this PPP, so SOPHIA may return 0 particles
    if (pnEnergy.size() == 0)
        return;

    double maxEnergy = *std::max_element(pnEnergy.begin(), pnEnergy.end());  // criterion for being declared primary
    for (int i = 0; i < pnEnergy.size(); ++i) {
		if (pnEnergy[i] == maxEnergy) {  // nucleon is primary particle
			if (A == 1) {
				// single interacting nucleon
				candidate->current.setEnergy(pnEnergy[i]);
				try
				{
					candidate->current.setId(sign * nucleusId(1, 14 - pnType[i]));
				}
				catch (std::runtime_error &e)
				{
					KISS_LOG_ERROR<< "Something went wrong in the PhotoPionProduction (primary particle, A==1)\n" << "Please report this error on https://github.com/CRPropa/CRPropa3/issues including your simulation setup and the following random seed:\n" << Random::instance().getSeed_base64();
					throw;
				}
			} else {
				// interacting nucleon is part of nucleus: it is emitted from the nucleus
				candidate->current.setEnergy(E - EpA);
				try
				{
					candidate->current.setId(sign * nucleusId(A - 1, Z - int(onProton)));
					candidate->addSecondary(sign * nucleusId(1, 14 - pnType[i]), pnEnergy[i], pos, tag);
				}
				catch (std::runtime_error &e)
				{
					KISS_LOG_ERROR<< "Something went wrong in the PhotoPionProduction (primary particle, A!=1)\n" << "Please report this error on https://github.com/CRPropa/CRPropa3/issues including your simulation setup and the following random seed:\n" << Random::instance().getSeed_base64();
					throw;
				}
			}
		} else {  // nucleon is secondary proton or neutron
			candidate->addSecondary(sign * nucleusId(1, 14 - pnType[i]), pnEnergy[i], pos, tag);
		}
	}
}

// double PhotoPionProduction::lossLength(int id, double gamma, double z) {
//     int A = massNumber(id);
//     int Z = chargeNumber(id);
//     int N = A - Z;

//     double lossRate = 0;
//     if (Z > 0)
//         lossRate += 1 / nucleonMFP(gamma, z, true) * nucleiModification(A, Z);
//     if (N > 0)
//         lossRate += 1 / nucleonMFP(gamma, z, false) * nucleiModification(A, N);

//     // approximate the relative energy loss
//     // - nucleons keep the fraction of mass to delta-resonance mass
//     // - nuclei lose the energy 1/A the interacting nucleon is carrying
//     double relativeEnergyLoss = (A == 1) ? 1 - 938. / 1232. : 1. / A;
//     lossRate *= relativeEnergyLoss;

//     // scaling factor: interaction rate --> energy loss rate
//     lossRate *= (1 + z);

//     return 1. / lossRate;
// }

} // namespace crpropa

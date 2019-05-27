#include "crpropa/module/HadronicInteraction.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"
#include "crpropa/Random.h"

#include <fstream>
#include <limits>
#include <cmath>
#include <stdexcept>

namespace crpropa {


HadronicInteraction::HadronicInteraction(double massDensity,
                                         bool electrons,
                                         bool photons,
                                         bool neutrinos,
                                         std::string tag) {
    setMassDensity(massDensity);
    this->spaceTimeGrid = ScalarGrid4d();
    this->spaceGrid = ScalarGrid();
    setHaveElectrons(electrons);
    setHavePhotons(photons);
    setHaveNeutrinos(neutrinos);
    this->tag = tag;
    setDescription("HadronicInteraction_isotropicConstant");
}

HadronicInteraction::HadronicInteraction(double massDensity,
                                         ScalarGrid4d spaceTimeGrid,
                                         bool electrons,
                                         bool photons,
                                         bool neutrinos,
                                         std::string tag) {
    setMassDensity(massDensity);
    this->spaceTimeGrid = spaceTimeGrid;
    this->spaceGrid = ScalarGrid();
    setHaveElectrons(electrons);
    setHavePhotons(photons);
    setHaveNeutrinos(neutrinos);
    this->tag = tag;
    setDescription("HadronicInteraction_spaceTimeDependent");
}

HadronicInteraction::HadronicInteraction(double massDensity,
                                         ScalarGrid spaceGrid,
                                         bool electrons,
                                         bool photons,
                                         bool neutrinos,
                                         std::string tag) {
    setMassDensity(massDensity);
    this->spaceTimeGrid = ScalarGrid4d();
    this->spaceGrid = spaceGrid;
    setHaveElectrons(electrons);
    setHavePhotons(photons);
    setHaveNeutrinos(neutrinos);
    this->tag = tag;
    setDescription("HadronicInteraction_spaceDependentConstant");
}

void HadronicInteraction::setMassDensity(double dens) {
    this->massDensity = dens;
}

void HadronicInteraction::setHaveElectrons(bool b) {
    this->haveElectrons = b;
}

void HadronicInteraction::setHavePhotons(bool b) {
    this->havePhotons = b;
}

void HadronicInteraction::setHaveNeutrinos(bool b) {
    this->haveNeutrinos = b;
}

Vector3d HadronicInteraction::getPosition(double height, double radius) const {
    Random &random = Random::instance();
    int i = 0;
    Vector3d pos(0, 0, 0);
    double phi = random.rand() * 2 * M_PI;
    int j = 0;
    do {
        double r = random.rand() * radius;
        double yr = random.rand();
        double Fr = exp(-r * r / (2 * 4.2 * 4.2 * kpc * kpc));
        if (yr < Fr) {
            pos = Vector3d(cos(phi) * r, sin(phi) * r, 0);
            j++;
        }
    } while (j == 0);
    do {
        double z = random.rand()* height;
        double yz = random.rand();
        double Fz = exp(-z / (10 * pc));
        if (yz < Fz) {
            double a = random.rand();
            if (a <= 0.5)
                z = -z;
            pos += Vector3d(0, 0, z);
            j++;
        }
    } while (j == 1);
    return pos;
}

double HadronicInteraction::distribution_e(double Eprimary, double x) const {
    /* Distribution function (energy) for electrons, electron neutrino and 
    (second) muon neutrinos based on Kelner 2006 - eqs. 62-65
    input: energy: primary's energy / x: Elepton/Epi | if second numu: x: Enumu/Eprimary
    */
    double L = log(Eprimary / TeV);
    double Be = 1 / (69.5 + 2.65 * L + 0.3 * L * L);
    double betae = 1 / pow((0.201 + 0.062 * L + 0.00042 * L * L), 0.25);
    double ke = (0.279 + 0.141 * L + 0.0172 * L * L) / (0.3 + (2.3 + L) * (2.3 + L));
    double F = Be * pow((1 + ke * pow(log(x), 2.)), 3.) / (x * (1 + 0.3 / pow(x, betae))) * (pow(-log(x), 5.));
    return F;
}

// Number of electrons, electron neutrinos and (second) muons neutrinos produced in a given interaction  based on Kelner 2006
int HadronicInteraction::numberOfElectrons(double Eprimary) const {
    const double xMax = 1.;
    const double xMin = 1. / 100000.;
    const double stepSize = 1. / 100000.;
    double x = xMin;
    double y = 0;
    double stepsDone = 0;
    do {
        y += distribution_e(Eprimary, x);
        x += stepSize;
        stepsDone++;
    } while (x < xMax);
    return round(y / stepsDone * (x - 1. / 1000.));
}

double HadronicInteraction::distribution_my1(double Eprimary, double x) const {
    /* Distribution function (energy) for (first) muon neutrino based on Kelner 2006 
    eqs. 66-69. input: energy: primary's energy / x: Enumu1/Eprimary*/
    double L = log(Eprimary / TeV);
    double Bm = 1.75 + 0.204 * L + 0.01 * pow(L, 2.);
    double betam = 1 / (1.67 + 0.111 * L + 0.0038 * pow(L, 2.));
    double km = 1.07 - 0.086 * L + 0.002 * pow(L, 2.);
    x /= 0.427;
    double aa = (1 - pow(x, betam)) / (1 + km * pow(x, betam) * (1 - pow(x, betam)));
    double A = Bm * log(x) / x * pow(aa, 4.);
    double B = 1 / log(x) - 4 * betam * pow(x, betam) / (1 - pow(x, betam))
             - 4 * km * betam * pow(x, betam) * (1 - 2 * pow(x, betam))
             / (1 + km * pow(x, betam) * (1 - pow(x, betam)));
    double F = A * B;
    return F;
}

// Number of (first) muon neutrinos produced in a given interaction  based on Kelner 2006
int HadronicInteraction::numberOfMuonNeutrinos(double Eprimary) const {
    const double xMax = 0.427;
    const double xMin = 1. / 100000.;
    const double stepSize = 1. / 100000.;
    double x = xMin;
    double y = 0.;
    int stepsDone = 0;
    do {
        y += distribution_my1(Eprimary, x);
        x += stepSize;
        stepsDone++;
    } while (x < xMax);
    return round(y / stepsDone * (x - 1. / 1000.));
}

double HadronicInteraction::distribution_gamma(double Eprimary, double x) const {
    /* Distribution function (energy) for gamma rays based on Kelner 2006 eqs. 58-61 
    energy: primary's energy / x: Egamma/Eprimary */
    double L = log(Eprimary / TeV);
    double Bg = 1.3 + 0.14 * L + 0.011 * L * L;
    double betag = 1 / (1.79 + 0.11 * L + 0.008 * L * L);
    double kg = 1 / (0.801 + 0.049 * L + 0.014 * L * L);
    double A = Bg * log(x) / x;
    double B = (1 - pow(x, betag)) / (1 + kg * pow(x, betag) * (1 - pow(x, betag)));
    double C = 1 / log(x) - 4 * betag * pow(x, betag) / (1 - pow(x, betag))
             - 4 * kg * betag * pow(x, betag) * (1 - 2 * pow(x, betag))
             / (1 + kg * pow(x, betag) * (1 - pow(x, betag)));
    double F = A * pow(B, 4.) * C;
    return F;
}

// Number of gamma rays produced in a given interaction  based on Kelner 2006
int HadronicInteraction::numberOfGammaRays(double Eprimary) const {
    const double xMin = 1. / 100000.;
    const double xMax = 1.;
    const double stepSize = 1. / 100000.;
    double x = xMin;
    double y = 0.;
    int stepsDone = 0;
    do {
        y += distribution_gamma(Eprimary, x);
        x += stepSize;
        stepsDone++;
    } while (x < xMax);
    return round(y / stepsDone * (x - 1. / 1000.));
}

// Energy distribution for lepton secondaries of pp interactions based on Carceller 2017
double HadronicInteraction::distribution_Carceller(double Eprimary, double x, double jcap, double a0, double b0) const {
    double a = a0 * (1 + 0.073 * log(Eprimary / PeV) + 0.0070 * log(Eprimary / PeV) * log(Eprimary / PeV));
    double b = b0 * (1 + 0.020 * log(Eprimary / PeV) + 0.0018 * log(Eprimary / PeV) * log(Eprimary / PeV));
    double A = a * pow((1 - jcap * x), 3.) / x;
    double B = exp(-b * pow(jcap * x, 0.43)) / pow(1 + pow(0.1 * GeV / (x * Eprimary), 0.5), 2.);
    double F = A * B;
    return F;
}

// Energy distribution for gamma photons based on Carceller 2017
double HadronicInteraction::distribution_Carceller_g(double Eprimary, double x, double jcap, double a0, double b0) const {
    double a = a0 * (1 + 0.073 * log(Eprimary / PeV) + 0.0070 * log(Eprimary / PeV) * log(Eprimary / PeV));
    double b = b0 * (1 + 0.020 * log(Eprimary / PeV) + 0.0018 * log(Eprimary / PeV) * log(Eprimary / PeV));
    double A = a * pow((1 - jcap * x), 3.) / x;
    double B = exp(-b * pow(jcap * x, 0.43)) / pow(1 + pow(0.2 * GeV / (x * Eprimary), 0.5), 2.);
    double F = A * B;
    return F;
}

// Cross Section of inelastic pp interaction based on Tan & Ng 1983 (Used in Galprop)
double HadronicInteraction::CrossSection_Galprop(double Eprimary) const {
    double cs_inel;
    double U = log(Eprimary / GeV * 1 / 200);
    if (U >= 0 and Eprimary >= 3 * GeV)
        cs_inel = (32.2 * (1 + 0.0273 * U)) * 1e-31 + 32.2 * 0.01 * pow(U, 2.) * 1e-31;
    if (U < 0 and Eprimary >= 3 * GeV)
        cs_inel = (32.2 * (1 + 0.0273 * U)) * 1e-31;
    if (Eprimary <= 0.3 * GeV)
        cs_inel = 0;
    return cs_inel;
}

// Cross Section of inelastic pp interaction based on Kelner 2006
double HadronicInteraction::CrossSection_Kelner(double Eprimary) const {
    double L = log(Eprimary / TeV);
    double A = 1 - pow(1.22 * 1e-3 * TeV / Eprimary, 4.);
    double cs_inel = (34.3 + 1.88 * L + 0.25 * L * L) * A * A * 1e-31;
    return cs_inel;
}

// Cross Section of inelastic pp interaction based on Carceller 2017
double HadronicInteraction::CrossSection_Carceller(double Eprimary) const {
    double cs_inel = 17.7 * pow(Eprimary / GeV, 0.082) * 1e-31;
    return cs_inel;
}

void HadronicInteraction::process(Candidate *candidate) const {
    // Interaction only for protons
    if (candidate->current.getId() != 1000010010)
        return;

    // Probability of interaction
    const double step = candidate->getCurrentStep();
    double Eprimary = candidate->current.getEnergy();
    const double cs_inel = CrossSection_Kelner(Eprimary);

    Vector3d pos = candidate->current.getPosition();
    const double time = candidate->getTrajectoryLength()/c_light;

    double dens = massDensity;
    const std::string description = getDescription();
    if (description == "HadronicInteraction_isotropicConstant") {
        // do nothing, just check for correct initialization
    } else if (description == "HadronicInteraction_spaceDependentConstant") {
        dens *= spaceGrid.interpolate(pos);
    } else if (description == "HadronicInteraction_spaceTimeDependent") {
        dens *= spaceTimeGrid.interpolate(pos, time);
    } else {
        throw std::runtime_error("HadronicInteraction: invalid description string");
    }
    const double p_pp = cs_inel * dens * step;

    // limit next step to mean free path
    const double limit = 1 / p_pp * 0.1;

    if (step > limit)
        candidate->limitNextStep(limit);

    // Interaction?
    Random &random = Random::instance();
    if (random.rand() > p_pp or Eprimary < 1 * GeV)
        return;

    pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());  // secondaries' pos

    /* Initialize current energies of secondaries */
    double Eout = 0;

    double Egamma = 0;  // gamma
    double EnuMu1 = 0;  // nu_mu1
    double Eelectron = 0;  // electron
    double EnuE = 0;  // nu_e
    double EnuMu2 = 0;  // nu_mu2
    double Etot = 0;

    /* Establish number of secondaries */
    const int Ngamma = numberOfGammaRays(Eprimary);  // number of gamma rays
    const int NnuMu1 = numberOfMuonNeutrinos(Eprimary);  // number of first myon neutrino
    const int Nelectron = numberOfElectrons(Eprimary);  // number of electron
    const int NnuE = Nelectron;  // number of electron neutrino
    const int NnuMu2 = Nelectron;  // number of second muon neutrino
    const int Ntotal = Nelectron + Ngamma + NnuMu1 + NnuE + NnuMu2;  // Total number of secondaries in the interaction

    /* Initialization for stopping criteria.
    Counter of interim produced particles: */
    int iTotal = 1;

    int iGamma = 1;
    int iNuMu1 = 1;
    int iElectron = 1;
    int iNuE = 1;
    int iNuMu2 = 1;

    const double threshold = 0.0001;
    int test;
    int end = -1;

    /* Jump to random starting point: */
    const double startPoint = random.rand();
    if (startPoint <= 0.2)
        goto label1;
    if (startPoint <= 0.4)
        goto label2;
    if (startPoint <= 0.6)
        goto label3;
    if (startPoint <= 0.8)
        goto label4;
    if (startPoint <= 1.0)
        goto label5;

    do {
        label1:

        // Gamma rays
        test = iGamma;
        // Check if all gamma rays were created
        if (iGamma <= Ngamma) {
            if (end == -1) {
                // pick gamma ray's energy
                do {
                    double x = threshold + random.rand() * (1 - Eout / Eprimary - threshold);
                    Eout = x * Eprimary;
                    double E = distribution_gamma(Eprimary, x);
                    double Emax = distribution_gamma(Eprimary, threshold);
                    double y = random.rand() * Emax;

                    if (y < E and (Etot + Eout) < Eprimary) {
                        if (havePhotons) {
                            if (1. / Eout != 0.)  // BUG: some photons are produced with infinite energy!
                        	    candidate->addSecondary(22, Eout, pos, tag);
                        }
                        Egamma += Eout;
                        Etot += Eout;
                        iTotal++;
                        iGamma++;
                        if (Etot / Eprimary >= (1 - threshold))
                            end = iTotal;
                    }
                } while (test == iGamma);
            } else {
                Eout = (Eprimary - Etot) / (Ntotal - end);
                if (havePhotons) {
                    if (1. / Eout != 0.)  // BUG: some photons are produced with infinite energy!
                        candidate->addSecondary(22, Eout, pos, tag);
                }
                Egamma += Eout;
                iTotal++;
                iGamma++;
            }
        }

        label2:

        // First myon neutrino 14
        test = iNuMu1;
        if (iNuMu1 <= NnuMu1) {
           if (end == -1) {
                do {
                    double x = threshold + random.rand() * (0.427 - threshold);
                    Eout = x * Eprimary;
                    double E = distribution_my1(Eprimary, x);
                    double Emax = distribution_my1(Eprimary, threshold);
                    double y = random.rand() * Emax;
                    if (y < E and (Etot + Eout) < Eprimary) {
                        if (haveNeutrinos)
                            candidate->addSecondary(14, Eout, pos, tag);
                        EnuMu1 += Eout;
                        Etot += Eout;
                        iTotal++;
                        iNuMu1++;
                        if (Etot / Eprimary >= (1 - threshold))
                            end = iTotal;
                    }
                } while (test == iNuMu1);
            } else {
                Eout = (Eprimary - Etot) / (Ntotal - end);
                if (haveNeutrinos)
                    candidate->addSecondary(14, Eout, pos, tag);
                EnuMu1 += Eout;
                iTotal++;
                iNuMu1++;
            }
        }

        label3:

        // Electron 11
        test = iElectron;
        if (iElectron <= Nelectron) {
            if (end == -1) {
                do {
                    double x = threshold + random.rand() * (1 - Eout / Eprimary - threshold);
                    Eout = x * Eprimary;
                    double E = distribution_e(Eprimary, x);
                    double Emax = distribution_e(Eprimary, threshold);
                    double y = random.rand() * Emax;
                    if (y < E and (Etot + Eout) < Eprimary) {
                        if (haveElectrons)
                            candidate->addSecondary(11, Eout, pos, tag);
                        Etot += Eout;
                        Eelectron += Eout;
                        iTotal++;
                        iElectron++;
                        if (Etot / Eprimary >= (1 - threshold))
                            end = iTotal;
                    }
                } while (test == iElectron);
            } else {
                Eout = (Eprimary - Etot) / (Ntotal - end);
                if (haveElectrons)
                    candidate->addSecondary(11, Eout, pos, tag);
                iTotal++;
                iElectron++;
                Eelectron += Eout;
            }
        }

        label4:

        // Electron neutrino 12
        test = iNuE;
        if (iNuE <= NnuE) {
            if (end == -1) {
                do {
                    double x = threshold + random.rand() * (1 - Eout / Eprimary - threshold);
                    Eout = x * Eprimary;
                    double E = distribution_e(Eprimary, x);
                    double Emax = distribution_e(Eprimary, threshold);
                    double y = random.rand() * Emax;
                    if (y < E and (Etot + Eout) < Eprimary) {
                        if (haveNeutrinos)
                            candidate->addSecondary(12, Eout, pos, tag);
                        EnuE += Eout;
                        Etot += Eout;
                        iTotal++;
                        iNuE++;
                        if (Etot / Eprimary >= (1 - threshold))
                            end = iTotal;
                    }
                } while (iNuE == test);
            } else {
                Eout = (Eprimary - Etot) / (Ntotal - end);
                if (haveNeutrinos)
                    candidate->addSecondary(12, Eout, pos, tag);
                iTotal++;
                iNuE++;
                EnuE += Eout;
            }
        }

        label5:

        // Second myon neutrino 14
        test = iNuMu2;
        if (iNuMu2 <= NnuMu2) {
            if (end == -1) {
                do {
                    double x = threshold + random.rand() * (1 - Eout / Eprimary - threshold);
                    Eout = x * Eprimary;
                    double E = distribution_e(Eprimary, x);
                    double Emax = distribution_e(Eprimary, threshold);
                    double y = random.rand() * Emax;
                    if (y < E and (Etot + Eout) < Eprimary) {
                        if (haveNeutrinos)
                            candidate->addSecondary(14, Eout, pos, tag);
                        EnuMu2 += Eout;
                        Etot += Eout;
                        iTotal++;
                        iNuMu2++;
                        if (Etot / Eprimary >= (1 - threshold))
                            end = iTotal;
                    }
                } while (iNuMu2 == test);
            } else {
                Eout = (Eprimary - Etot) / (Ntotal - end);
                if (haveNeutrinos)
                    candidate->addSecondary(14, Eout, pos, tag);
                iTotal++;
                iNuMu2++;
                EnuMu2 += Eout;
            }
        }
    } while (iTotal <= Ntotal);

    if (end != -1)
        std::cout << end << " end != -1" << std::endl;

    // Reduce primary's energy
    Eprimary -= (EnuE + EnuMu2 + Eelectron + EnuMu1 + Egamma);
    if (Eprimary <= 0.) {
        std::cout << "warning: Eprimary = " << Eprimary / GeV <<  " GeV" << std::endl;
        candidate->setActive(false);
    }
    candidate->current.setEnergy(Eprimary);
    return;
}

} // namespace CRPropa

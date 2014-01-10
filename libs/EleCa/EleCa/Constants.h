#ifndef ELECA_CONSTANTS_H_
#define ELECA_CONSTANTS_H_

#include <cmath>

namespace eleca {

static const double cPI = M_PI;

static const int POINTS_VERY_FEW = 100;

static const double M2MPC = 3.240779e-17 / 1.0e6;
static const double KM2MPC = 3.240779e-20;
static const double S2YR = 3.168808781403e-08;
static const double H0 = 70.4;
static const double H0y = H0 * (double) KM2MPC / S2YR;   // H0 per years
static const double OC = 0.227;  // Dark matter density
static const double OB = 0.0456; // Baryon density
static const double OM = OB + OC;  // Matter density
static const double OL = 0.728;  // Dark energy density

static const double K_CBR = 1.318684673251832e+13; // =1/pi**2/hbarc**3 [eV^-3 cm^-3]
static const double eV2J = 1.602176487e-19; // from eV to J
static const double ElectronMass = 0.510998918e6; // [eV/c^2]
static const double K_boltz = 8.617342294984e-5;  // [eV/K ] Boltzman constant
static const double C_speed = 299792458;          // [m/s] speed of light
static const double SigmaThompson = 6.6524e-25;

static const double T_CMB = 2.725; // [K] // evolution 2.725*(1-z)  1012.3164 
static const double T_COB = 5270; // [K]  // Visible [380 - 760] nm. Scelgo 550
static const double T_CIB = 1.45e+02; // [k] Middle IR 5 to (25-40) µm according to Nasa. scelgo 20e-6 m
static const double T_CRB = 3e-03; // [k] ~ cm - 10m.  scelgo ~1 m

static const double CMB_en = K_boltz * T_CMB; //2.348e-4;             // [eV]
static const double CRB_en = K_boltz * T_CRB;
static const double COB_en = K_boltz * T_COB;
static const double CIB_en = K_boltz * T_CIB;
static const double CIOB_en = CIB_en + COB_en;    // [eV]

static const double h_Planck = 4.135667e-15; // [eV s]// 
static const double hcut_Planck = h_Planck / 2 / cPI; // [eV s] hcut = h/2Pi [Js]
static const double LambdaCompton = hcut_Planck / (ElectronMass / C_speed);

static const double eps_ph_inf_urb = 4.1e-12;   // [eV]
static const double eps_ph_inf_cmb = 0.825e-6;   // [eV]
static const double eps_ph_inf_cib = 2e-3;   // [eV]
static const double eps_ph_inf_cob = 5e-2;   // [eV]
static const double eps_ph_inf_ciob = 2e-3;   // [eV]

static const double eps_ph_sup_urb = eps_ph_inf_cmb;   //4e-5;   // [eV]
static const double eps_ph_sup_cmb = eps_ph_inf_cob;   // [eV]
static const double eps_ph_sup_cob = 9.9;   // [eV]
static const double eps_ph_sup_cib = 0.8;   // [eV]
static const double eps_ph_sup_ciob = 9.9;   // [eV]

static const double eps_ph_sup_global = eps_ph_sup_cob;   // [eV] *global
static const double eps_ph_inf_global = eps_ph_inf_urb;   // [eV] *global

static const int NsecG = 0;

static const double z0ph = 0;
static const double E0ph = 0;

static const int particle_type = 0;
static const double EnergyCM = 0;
static const double GammaEnergy = 0;
static const double BackGamma = 0;
static const double PPxsection = 0;
static const double DPPxsection = 0;
static const double TPPxsection = 0;
static const double ICSxsection = 0;
static const double PPlength = 0;
static const double DPPlength = 0;
static const double TPPlength = 0;
static const double ICSlength = 0;
static const double n_eps2 = 0;
static const double eps2 = 0;

} // namespace eleca

#endif // ELECA_CONSTANTS_H_

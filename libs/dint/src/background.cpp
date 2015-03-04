
#include "dint/background.h"
#include <cmath>

// integer pow implementation as template that is evaluated at compile time
template <unsigned int exponent>
inline double pow_integer(double base)
{
  return pow_integer<exponent >> 1>(base*base) * (((exponent & 1) > 0) ? base : 1);
}

template <>
inline double pow_integer<0>(double base)
{
  return 1;
}



void LoadPhotonBackground(const double redshift, 
                          dCVector* pBgEnergy, dCVector* pBgEnergyWidth, 
                          dCVector* pBgPhotonDensity,
			  const int aIRFlag, const double aZmax_IR, const int aRadioFlag,
			  const double aH0, const double aOmegaM, const double aOmegaLambda) {
  if ((pBgEnergy->dimension != pBgEnergyWidth->dimension) ||
      (pBgEnergyWidth->dimension != pBgPhotonDensity->dimension))
    Error("LoadPhotonBackground: inconsistent dimensions", PROGRAM_ERROR);
  
  LoadCMB(redshift, pBgEnergy, pBgEnergyWidth, pBgPhotonDensity);
  LoadIR(redshift, pBgEnergy, pBgEnergyWidth, pBgPhotonDensity, aIRFlag, aZmax_IR);
  LoadRadio(redshift, pBgEnergy, pBgEnergyWidth, pBgPhotonDensity, aRadioFlag,
    	    aH0, aOmegaM, aOmegaLambda);
  
  for (int i=0; i<pBgEnergy->dimension; i++)
    (pBgPhotonDensity->vector)[i] *= VOLUME_UNIT; // normalization

#ifdef PRINT_PHOTON_BACKGROUND  // Output the density and exit (July 2005)
  cout << "z=" << redshift << endl;
  for (int i=0; i<pBgEnergy->dimension; i++)
    cout << (pBgEnergy->vector)[i] << " " << (pBgPhotonDensity->vector)[i] << endl;
  exit(-1) ;
#endif

}


void LoadCMB(const double redshift, const dCVector* pBgEnergy, 
             const dCVector* pBgEnergyWidth, dCVector* pBgPhotonDensity) {
  double numberDensity = 0.;    // to test CMB was computed correctly
  const double PRESENT_TEMPERATURE = 2.73; // present CMB temperature (K)
  double exponent;
  double tempNumber;
  double temperature;
  
  temperature = 1.6863e-10*PRESENT_TEMPERATURE*(1. + redshift); 
  //temperature at redshift (K)
  
  for (int i=0; i<pBgEnergy->dimension; i++) {
    exponent = (pBgEnergy->vector)[i]/temperature;
    if (exponent < 550.) {
      tempNumber = (pBgEnergy->vector)[i]*(pBgEnergy->vector)[i]/
	(exp(exponent) - 1.); // thermal distribution
    } else tempNumber = 0.;
    (pBgPhotonDensity->vector)[i] = tempNumber*
      (pBgEnergyWidth->vector)[i]*1.75586e30;
    numberDensity += (pBgPhotonDensity->vector)[i];
  }

#ifdef TEST_CMB
  printf("Computed CMB number density: %g\n", numberDensity);
  numberDensity = 20.2417*PRESENT_TEMPERATURE*PRESENT_TEMPERATURE*
    PRESENT_TEMPERATURE*(1. + redshift)*(1. + redshift)*
    (1. + redshift);
  printf("Predicted CMB number density: %g\n", numberDensity);
#endif
  
}


double CIOBR(double eps) {
	// parametrization for infrared/optical by
	// Hopkins, A. M. & Beacom, J. F. 2006, ApJ, 651, 142
	// See Model D Finke et al, arXiv:0905.1115v2
	const double eps_ph_sup_ciob = 9.9;   // [eV]
	const double eps_ph_inf_ciob = 2e-3;   // [eV]

	double tmp = 0;

	if (eps > eps_ph_inf_ciob && eps < eps_ph_sup_ciob) {
		double x = log(eps);
		tmp = -5.32524895349885 - 0.0741140642891119 * x
				- 0.252586527659431 * pow_integer<2>(x)
				+ 0.234971297531891 * pow_integer<3>(x)
				- 0.217014471117521 * pow_integer<4>(x)
				- 0.364936722063572 * pow_integer<5>(x)
				+ 0.0880702191711222 * pow_integer<6>(x)
				+ 0.221947767409286 * pow_integer<7>(x)
				+ 0.0445499623085708 * pow_integer<8>(x)
				- 0.0517435600939147 * pow_integer<9>(x)
				- 0.0295646851279071 * pow_integer<10>(x)
				- 0.00011943632049331 * pow_integer<11>(x)
				+ 0.00461621589174355 * pow_integer<12>(x)
				+ 0.00150906100702171 * pow_integer<13>(x)
				+ 1.91459088023263e-05 * pow_integer<14>(x)
				- 0.000110272619218937 * pow_integer<15>(x)
				- 3.45221358079085e-05 * pow_integer<16>(x)
				- 5.42000122025042e-06 * pow_integer<17>(x)
				- 4.90862622314226e-07 * pow_integer<18>(x)
				- 2.45145316799091e-08 * pow_integer<19>(x)
				- 5.25792204884819e-10 * pow_integer<20>(x);
		tmp = 0.4 * (double) exp(tmp) / eps / eps;
	} else {
		tmp = 0;
	}

	if (std::isnan(tmp))
		tmp = 0;

	return tmp;
}

double CIB_Evolution_Fast(double z) {
	// Function for the CIB fast evolution.
	// Stecker, Malkan, Scully (2006) arXiv:astro-ph/0510449v4

	double tmp = 0;
	double m = 4.;
	double z_flat = 1.;

	if (z <= z_flat)
		tmp = pow(1. + z, m);
	if (z_flat < z && z < 6)
		tmp = pow(1. + z_flat, m);
	if (z > 6)
		tmp = 0;

	return tmp;
}

double CIOB_Evolution(double z) {
	return CIB_Evolution_Fast(z);
}

double ElecaIOBR(double E, double z)
{
	double eps = E * ELECTRON_MASS;
	return CIOBR(eps) * CIOB_Evolution(z) * ELECTRON_MASS;
}


void LoadIR(const double redshift, const dCVector* pBgEnergy,
            const dCVector* pBgEnergyWidth, dCVector* pBgPhotonDensity,
	    const int aIRFlag, const double aZmax_IR) {

  double IRStartRedshift = aZmax_IR ;  // Primack CDM
  //  const double IRStartRedshift = 4.3; // Primack CHDM
  // parameters from Franceschini, Yoshi, and Takahara burst of both components
  const double deltaO = 7.;
  const double deltaD = 7.;
  double (*IRFunction)(const double zTarget, const double zObserve,
		       const double energy0, const double deltaO,
		       const double deltaD);
  if (aIRFlag == 0) IRFunction = HighIR ;
  else if (aIRFlag == 1) IRFunction = LowIR ;
  else if (aIRFlag == 2) IRFunction = NULL ;
  else if (aIRFlag == 4) IRFunction = NULL;
  else Error("LoadIR : Uncorrect IR flag",IO_ERROR) ;

	// eleca IR model
	if (aIRFlag == 4)
	{
    for (int i=0; i<pBgEnergy->dimension; i++) 
		{
		
				(pBgPhotonDensity->vector)[i] += ElecaIOBR((pBgEnergy->vector)[i], redshift)*(pBgEnergyWidth->vector)[i];
		}
	
	}
	else if (redshift <= IRStartRedshift) 
	{
    double z[GAULEG_POINTS];
    double w[GAULEG_POINTS];
    double observeRedshift;
    for (int i=0; i<pBgEnergy->dimension; i++) 
		{
      if (aIRFlag == 2) {
	(pBgPhotonDensity->vector)[i] += 
	  IR2(redshift,(pBgEnergy->vector)[i])*(pBgEnergyWidth->vector)[i]*ELECTRON_MASS;
      } else {
				double temp = 0;
				Gauleg(redshift, IRStartRedshift, z, w, GAULEG_POINTS); //find integration points
				observeRedshift = redshift;
				for (int j=0; j<GAULEG_POINTS; j++)
				{
					temp += IRFunction(z[j], observeRedshift, 
						(pBgEnergy->vector)[i], deltaO, deltaD)*w[j];
				}
				(pBgPhotonDensity->vector)[i] += temp*(pBgEnergyWidth->vector)[i];
      }
    }
  }

}


double IR2(const double redshift, const double BgEnergy) { 

  // Primack et al. (1999) CIB parameters
 
  const double flux_conversion = 2.9979e10/4./PI*1.602e-19*1.e9*1.e4;
  // xData is log10(wavelength/micrometers)
  // yData is log10(\nu I_\nu/[nW m^-2 sr^-1])

  // Redshift effect : we assume the following evolution of photon number density :
  //  n(e,z) = (1+z)**2 * n_0(e/(1+z))
  // It's the same as for CMB (pysically no injection of photons) and therefore the total 
  // number density evolves as (1+z)**3

  // Changed to Kneiske Background for consistency with CRPropa interactions, JK
  // Kneiske et al. Astron.\ Astrophys.\  {\bf 413} (2004) 807 
  const double xData[15] = {-1.00000, -0.750000, -0.500000, -0.250000, 0.00000, 0.250000, 0.500000, 0.750000,
			    1.00000, 1.25000, 1.50000, 1.75000, 2.00000, 2.25000, 2.50000};
  const double yData[15] = { -0.214401, 0.349313, 0.720354, 0.890389, 1.16042, 1.24692, 1.06525, 0.668659, 0.536312, 0.595859, 0.457456,
			     0.623521, 1.20208, 1.33657, 1.04461};

//   const double xData[15] = {-1., -0.75, -0.5, -0.25, 0., 0.25, 0.5,
// 			    0.75, 1., 1.25, 1.5, 1.75, 2., 2.25, 2.5};
//   const double yData[15] = {0.8, 1.1, 1.15, 1.2, 1.3, 1.2, 1.05, 0.7,
// 			    0.4, 0.3, 0.5, 0.8, 1.1, 1.3, 1.};

  double result;
  double lLocalEnergy = BgEnergy / (1.+redshift) ; // Redshift effect
  double x = C/(lLocalEnergy*ELECTRON_MASS*2.42e14)*1.e4; // convert energy into microns
  if (x > 500. || log10(x) <= xData[0]) {
    result = 0.;
  } else if (log10(x) >= xData[14]) {
    result = (yData[14] - yData[13])/(xData[14] - xData[13])*
      (log10(x) - xData[13]) + yData[13];
    result = pow(10., result);
  } else {
    int index = 1;
    while (xData[index] < log10(x)) index++;
    result = (yData[index] - yData[index-1])/
      (xData[index] - xData[index-1])*(log10(x) - xData[index-1]) + 
      yData[index-1];
    result = pow(10., result);
  }
  result *= pow(1.+redshift,2.)/pow(lLocalEnergy*ELECTRON_MASS, 2.)/
    flux_conversion; // Redshift effect

  return result;
}


double HighIR(const double zTarget, const double zObserve, 
	      const double energy0, const double deltaO,
	      const double deltaD) {
  const double normalization = 2./125.*1.e-5;
  // high norm: this is the only difference from LowIR

  double coefficient = 1./pow(1. + zTarget, 4.5)/H_0*4.*PI*
    (1. + zObserve)*(1. + zObserve);
  double energyTarget = energy0*(1. + zTarget)/(1. + zObserve);
  
  double function = OpticalIR(energyTarget)*pow(1. + zTarget, deltaO)/30. +
    DustIR(energyTarget)/300.*pow(1. + zTarget, deltaD)/2.;
  // for no evolution for both (i.e. initial burst)
  
  return (normalization*coefficient*function) ;
}


double LowIR(const double zTarget, const double zObserve, 
	     const double energy0, const double deltaO,
	     const double deltaD) {
  const double normalization = 0.4/125.*1.e-5;
  // low norm: this is the only difference from HighIR

  double coefficient = 1./pow(1. + zTarget, 4.5)/H_0*4.*PI*
    (1. + zObserve)*(1. + zObserve);
  double energyTarget = energy0*(1. + zTarget)/(1. + zObserve);
  
  double function = OpticalIR(energyTarget)*pow(1. + zTarget, deltaO)/30. +
    DustIR(energyTarget)/300.*pow(1. + zTarget, deltaD)/2.;
  // for no evolution for both (i.e. initial burst)
 
  return (normalization*coefficient*function) ;
}


double OpticalIR(const double energy) {
  const double xData[9] = {-1.0402, -0.9, -0.7, -0.63, -0.4, -0.1, 0.1, 0.6,
			   1.};
  const double yData[9] = {-30., -1.45, -1.65, -1.74, -0.33, 0.1, 0.2, -0.6, 
			   -2.6};
  double result;
  double x = C/(energy*ELECTRON_MASS*2.42e14)*1.e4; // convert energy into microns
  // cut off far IR excess completely

  if (x > 100. || log10(x) <= xData[0]) {
    result = 0.;
  } else if (log10(x) >= xData[8]) {
    result = (yData[8] - yData[7])/(xData[8] - xData[7])*
      (log10(x) - xData[7]) + yData[7];
    result = pow(10., result)/energy/energy;
  } else {
    int index = 1 ;
    while (xData[index] < log10(x)) index++;
    result = (yData[index] - yData[index-1])/
      (xData[index] - xData[index-1])*(log10(x) - xData[index-1]) + 
      yData[index-1];
    result = pow(10., result)/energy/energy;
  }
  
  return result;
}


double DustIR(const double energy) {
  const double xData[6] = {-0.18, 0.54, 1.06, 1.4, 1.8, 2.9};
  const double yData[6] = {-3.0, 0.77, 0.6, 1.37, 1.5, -3.0};
  double result;
  double x = C/(energy*ELECTRON_MASS*2.42e14)*1.e4;

  if (log10(x) <= xData[0]) result = 0.;
  else if (log10(x) >= xData[5]) {
    result = (yData[5] - yData[4])/(xData[5] - xData[4])*
      (log10(x) - xData[4]) + yData[4];
    result = pow(10., result)/energy/energy;
  } else {
    int index = 1 ;
    while (xData[index] < log10(x)) index++;
    result = (yData[index] - yData[index-1])/
      (xData[index] - xData[index-1])*(log10(x) - xData[index-1]) + 
      yData[index-1];
    result = pow(10., result)/energy/energy;
  }

  return result;
}


// Universal Radio Background from Protheroe, Bierman 1996. 
double URB(double eps) {
	const double eps_ph_inf_urb = 4.1e-12;   // [eV]
	const double eps_ph_inf_cmb = 0.825e-6;   // [eV]
	const double eps_ph_sup_urb = 4E-5;   //4e-5;   // [eV]
	if (eps < eps_ph_inf_urb || eps > eps_ph_sup_urb)
		return 0;

	const double h_Planck = 4.135667e-15; // [eV s]// 
	const double C_speed = 299792458;          // [m/s] speed of light
	const double eV2J = 1.602176487e-19; // from eV to J

	double v = eps / h_Planck;
	double x = log10(v / 1e9);

	double p0 = -2.23791e+01;
	double p1 = -2.59696e-01;
	double p2 = 3.51067e-01;
	double p3 = -6.80104e-02;
	double p4 = 5.82003e-01;
	double p5 = -2.00075e+00;
	double p6 = -1.35259e+00;
	double p7 = -7.12112e-01;  //xbreak

	double intensity = 0;
	if (x > p7)
		intensity = p0 + p1 * x + p3 * x * x * x / (exp(p4 * x) - 1) + p6 + p5 * x;
	else
		intensity = p0 + p1 * x + p2 * x * x
				+ p3 * x * x * x / (exp(p4 * x) - 1);
	intensity = pow(10, intensity);
	double n_eps = 0;
	n_eps = 4 * M_PI / (h_Planck * C_speed) * (intensity / eps);
	return n_eps / eV2J / 1.0e6;
}
double URB_Evolution(double z) {
	//from Protheroe - Bierman astro-ph:9605119
	if (z < 0.8)
		return pow_integer<4>(1. + z);
	return pow_integer<4>(1 + 0.8);   // z>= z0
}

/// Provide Radio Background from Eleca for Dint
double ElecaRadio(const double zTarget, const double zObserve, 
		const double energy0, const double aH0,
		const double aOmegaM, const double aOmegaLambda) {
	
	// Convert DINT unit to eps
	const double eps = energy0 * ELECTRON_MASS;
	return URB(eps) * URB_Evolution(zTarget) *  ELECTRON_MASS;
}


void LoadRadio(const double redshift, const dCVector* pBgEnergy,
               const dCVector* pBgEnergyWidth, dCVector* pBgPhotonDensity,
	       const int aRadioFlag, const double aH0, const double aOmegaM,
	       const double aOmegaLambda) {
  const double radioStartRedshift = 5.;
  // although the source distribution does not vanish at high redshift,
  //   one might as well cut it off effectively 
  //   (5 is 6 sigma away from peak)
  double (*RadioFunction)(const double zObserve, const double zTarget, 
			  const double energy0, const double H,
			  const double OM, const double OL);
  const double B = 2.58734e-28;
	if (aRadioFlag == 4)
	{
		for (int i=0; i<pBgEnergy->dimension; i++) {
					(pBgPhotonDensity->vector)[i] += 
						ElecaRadio(redshift, redshift, 
											(pBgEnergy->vector)[i], aH0, aOmegaM, aOmegaLambda) * (pBgEnergyWidth->vector)[i];
      }
	}
	else if (aRadioFlag != 3) {
    // aRadioFlag == 3 : Null radio, so do nothing
    if (aRadioFlag == 0) {
      RadioFunction = HighRadio ;
    } else if (aRadioFlag == 1) {
      RadioFunction = MedRadio ;
    } else if (aRadioFlag == 2) {
      RadioFunction = ObsRadio ;
    } else Error("LoadRadio : Uncorrect Radio flag",IO_ERROR) ;
    
    if (redshift <= radioStartRedshift  ) {
      double z[GAULEG_POINTS];
      double w[GAULEG_POINTS];
      double cutoffFactor; // cutoffFactor removed from function
      
      for (int i=0; i<pBgEnergy->dimension; i++) {
				double temp = 0 ;
				Gauleg(redshift, radioStartRedshift, z, w, GAULEG_POINTS);
				double observeRedshift = redshift;
				for (int j=0; j<GAULEG_POINTS; j++){
					temp += RadioFunction(z[j], observeRedshift, 
											(pBgEnergy->vector)[i], aH0, aOmegaM, aOmegaLambda)*w[j];
				}
					if (RadioFunction == ObsRadio)
							cutoffFactor = exp(-B/(pBgEnergy->vector)[i]/ (pBgEnergy->vector)[i]);
							// this is a cutoff of the spectrum due to ISM absorption
					else cutoffFactor = 1.;
					(pBgPhotonDensity->vector)[i] += cutoffFactor*temp * (pBgEnergyWidth->vector)[i];
      }
    }
  }
}


double HighRadio(const double zTarget, const double zObserve, 
		 const double energy0, const double aH0,
		 const double aOmegaM, const double aOmegaLambda) {
  double B0 = 25.26*log(10.);
  double B1 = 1.18*log(10.);
  double B2 = 0.28*log(10.);
  double H_z,coefficient,function,energyTarget;

  H_z = aH0*sqrt(aOmegaM*pow(1.+zTarget,3.)+aOmegaLambda);
  energyTarget = energy0*(1. + zTarget)/(1. + zObserve);
  coefficient = 4.*PI/H_z*DISTANCE_UNIT*1.e6/1.e5*4.2458e-56*
    pow(1. + zObserve, 2.)/(1.40576e-28 + 1.81075e-45/energyTarget + 
			    6.38227e-8*pow(energyTarget, 1.3) + pow(energyTarget, 1.8));
  function = exp(B0 + B1*zTarget - B2*zTarget*zTarget);
  
  return (1.e0*coefficient*function);
}


double MedRadio(const double zTarget, const double zObserve, 
		const double energy0, const double aH0,
		const double aOmegaM, const double aOmegaLambda) {
  
  return (HighRadio(zTarget,zObserve,energy0,aH0,aOmegaM,aOmegaLambda)/pow(10.,0.7));
}


double ObsRadio(const double zTarget, const double zObserve, 
		const double energy0, const double aH0,
		const double aOmegaM, const double aOmegaLambda) {
  double B0 = 25.26*log(10.);
  double B1 = 1.18*log(10.);
  double B2 = 0.28*log(10.);
  double H_z,coefficient,function,energyTarget;

  H_z = aH0*sqrt(aOmegaM*pow(1.+zTarget,3.)+aOmegaLambda);
  energyTarget = energy0*(1. + zTarget)/(1. + zObserve);
  coefficient = 4.*PI/H_z*DISTANCE_UNIT*1.e6/1.e5*4.2458e-56*
    pow(1. + zObserve, 2.)*pow(energyTarget, -1.8);
  function = exp(B0 + B1*zTarget - B2*zTarget*zTarget);
  
  return (coefficient*function);
}




/*
// Routine to "dump" the spectrum.
// Not needed
void DumpBgSpectrum(const dCVector* pBgEnergy, const dCVector* pBgEnergyWidth,
		    const dCVector* pBgPhotonDensity, const char* filename)
{
    int i;
    int j;
    FILE* dumpFile;
    char f1[80] = "datafiles/";

    if ((pBgEnergy->dimension != pBgEnergyWidth->dimension) ||
	(pBgEnergyWidth->dimension != pBgPhotonDensity->dimension))
    {
	Error("DumpBgSpectrum: inconsistent dimensions", PROGRAM_ERROR);
    }

    strncat(f1, filename, 79 - strlen(filename));
    dumpFile = SafeFOpen(f1, "w");
// this is to send the dump file to a different directory (datafiles)
//     by Guenter (7/20/1998) 

    for (i = 0; i < pBgEnergy->dimension; i++)
	{
	    fprintf(dumpFile, "%15.4E %15.4E\n", 
		    ELECTRON_MASS*(pBgEnergy->vector)[i],
		    (pBgPhotonDensity->vector)[i]/ELECTRON_MASS/
		    (pBgEnergyWidth->vector)[i]/VOLUME_UNIT);
// proper unit conversion for energy 
	}
    fclose(dumpFile);
}
*/

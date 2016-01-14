#ifndef ELECA_ENERGY_LOSS_H_
#define ELECA_ENERGY_LOSS_H_

#include "EleCa/EnergyLoss.h"
#include "EleCa/Process.h"
#include "EleCa/Common.h"

#include <iostream>
#include <stdexcept>
#include <cstdlib>
#include <cfloat>

namespace eleca {

static const int MC_SAMPLING = 1000;

double dTdZ(double z) {
	// See High Energy Cosmic Rays, Todor Stanev, Pag. 232 (2009)
	return -1.
			/ ((1 + z) * H0y
					* sqrt(
							pow_integer<3>(1. + z) * OM + OL
									+ (1 - OM - OL) * pow_integer<2>(1. + z)));
}

double betaRsh(double z) {
	// Energy loss term due to cosmological redshift
	return H0y
			* sqrt(pow_integer<3>(1. + z) * OM + OL + (1 - OM - OL) * pow_integer<2>(1. + z));
}

double fLossAdiabatic(double E, double z) {
	return -dTdZ(z) * betaRsh(z) * E;
}

double AdiabaticELoss(double z0, double z, double E0) {
	return E0 * (double) (1. + z) / (1. + z0);
}

// ##################################################
// # Synchrotron rate of Energy loss averaged over angles
// ##################################################
// for an ensemble of electrons that are scattered randomly in all directions:
double MeanRateSynchrotronLoss(double E, double B) {
	double dEdt = 0;
	if (B > 0)
		dEdt = 3.79e-6 * pow_integer<2>(E / 1e9 * B) * 1e9 / E;

	return dEdt;
}

// ####################################
// # Synchrotron Energy loss
// ####################################

double ESynchrotronAdiabaticLoss(double z, double E, double B) {
	double dEdt = MeanRateSynchrotronLoss(E, B);

	return E * (-dTdZ(z) * (dEdt + betaRsh(z)));
}

static const int RK_ORDER = 6;
static double gRKa[RK_ORDER + 1];
static double gRKc[RK_ORDER + 1];
static double gRKcs[RK_ORDER + 1];
static double gRKb[RK_ORDER + 1][RK_ORDER];
static bool gRKInitialized = false;

void InitRK() {
	if (gRKInitialized)
		return;

	// Current Runge-Kutta method for solving ODE
	gRKa[0] = 0;
	gRKa[1] = 0;
	gRKa[2] = 1. / 5.;
	gRKa[3] = 3. / 10.;
	gRKa[4] = 3. / 5.;
	gRKa[5] = 1.;
	gRKa[6] = 7. / 8.;

	for (int i = 0; i < RK_ORDER + 1; i++) {
		for (int j = 0; j < RK_ORDER; j++)
    {
			gRKb[i][j] = 0.;
    }
	}

	gRKb[2][1] = 1. / 5.;
	gRKb[3][1] = 3. / 40.;
	gRKb[3][2] = 9. / 40.;
	gRKb[4][1] = 3. / 10.;
	gRKb[4][2] = -9. / 10.;
	gRKb[4][3] = 6. / 5.;
	gRKb[5][1] = -11. / 54.;
	gRKb[5][2] = 5. / 2.;
	gRKb[5][3] = -70. / 27.;
	gRKb[5][4] = 35. / 27.;
	gRKb[6][1] = 1631. / 55296.;
	gRKb[6][2] = 175. / 512.;
	gRKb[6][3] = 575. / 13824.;
	gRKb[6][4] = 44275. / 110592.;
	gRKb[6][5] = 253. / 4096.;

	gRKc[0] = 0.;
	gRKc[1] = 37. / 378.;
	gRKc[2] = 0.;
	gRKc[3] = 250. / 621.;
	gRKc[4] = 125. / 594.;
	gRKc[5] = 0.;
	gRKc[6] = 512. / 1771.;

	gRKcs[0] = 0.;
	gRKcs[1] = 2825. / 27648.;
	gRKcs[2] = 0.;
	gRKcs[3] = 18575. / 48384.;
	gRKcs[4] = 13525. / 55296.;
	gRKcs[5] = 277. / 14336.;
	gRKcs[6] = 1. / 4.;

	gRKInitialized = true;
}

//===================================

double EnergyLoss1D(double Energy, double z0, double zfin, double B) {

	double zStep = 2.5e-5;

#pragma omp critical
	{
		InitRK();
	}
	double k1, k2, k3, k4, k5, k6;

	bool FLAG_PROPAG = 1;

	while (z0 > zfin && FLAG_PROPAG) {

		k1 = -zStep * ESynchrotronAdiabaticLoss(z0 - zStep, Energy, B);
		k2 = -zStep
				* ESynchrotronAdiabaticLoss(z0 - zStep * gRKa[2],
						Energy + k1 * gRKb[2][1], B);

		k3 = -zStep
				* ESynchrotronAdiabaticLoss(z0 - zStep * gRKa[3],
						Energy + k1 * gRKb[3][1] + k2 * gRKb[3][2], B);
		k4 = -zStep
				* ESynchrotronAdiabaticLoss(z0 - zStep * gRKa[4],
						Energy + k1 * gRKb[4][1] + k2 * gRKb[4][2]
								+ k3 * gRKb[4][3], B);
		k5 = -zStep
				* ESynchrotronAdiabaticLoss(z0 - zStep * gRKa[5],
						Energy + k1 * gRKb[5][1] + k2 * gRKb[5][2]
								+ k3 * gRKb[5][3] + k4 * gRKb[5][4], B);
		k6 = -zStep
				* ESynchrotronAdiabaticLoss(z0 - zStep * gRKa[6],
						Energy + k1 * gRKb[6][2] + k2 * gRKb[6][2]
								+ k3 * gRKb[6][3] + k4 * gRKb[6][4]
								+ k5 * gRKb[6][5], B);

		Energy = Energy
				+ (k1 * gRKc[1] + k2 * gRKc[2] + k3 * gRKc[3] + k4 * gRKc[4]
						+ k5 * gRKc[5] + k6 * gRKc[6]);

		z0 -= zStep;

		if (fabs(z0) < 1e-8 || z0 < 0) {
			z0 = 0.;
			FLAG_PROPAG = 0;
		}
		if (fabs(z0 - zfin) < 1e-8 || z0 < zfin)
			z0 = zfin;

		if (Energy < 1e9) {
			FLAG_PROPAG = 0;
		}

		if (std::isnan(Energy))
			return 0;
	}

	return Energy;
}

double dSigmadE_ICS(double Ee, double Eer, double s, double theta) {
	/*!
	 Differential cross-section for inverse Compton scattering. from lee, eq. 23
	 */

	double beta = (s - ElectronMass * ElectronMass)
			/ (s + ElectronMass * ElectronMass);
	// boundaries rewritten to avoid error due to numerical uncertainties
	if ((1 - Eer / Ee) / (Eer / Ee +1) - beta > DBL_EPSILON || Eer / Ee > 1) 
	{
		std::cerr << "ERROR, Energy outside limits for ICS [Lee96]! " << std::endl;
		std::cerr << "       Eer = " << Eer << " Ee = " << Ee << "  Eer/Ee = " <<
			Eer / Ee << "   (1 - beta) / (1 + beta) = " << (1 - beta) / (1 + beta) <<
			" beta = " << beta << std::endl;
		return 0.;
	}
	else
	{
		double q = ((1 - beta) / beta) * (1 - Ee / Eer);
		double A = Eer / Ee + Ee / Eer;
		double k = (3.0 / 8.0) * (SigmaThompson * ElectronMass * ElectronMass)
				/ (s * Ee);
		double dsigmadE = k * ((1 + beta) / beta) * (A + 2 * q + q * q);

		return dsigmadE;
	}
}

  double dSigmadE_PP(double Ee, double E0, double eps, double theta, double s) {
	/*!
	 Differential cross-section for pair production.
	 */
    //	double s = ElectronMass * ElectronMass + 2 * eps * E0 * (1 - cos(theta));
	double beta = sqrt(1 - 4 * ElectronMass * ElectronMass / s);

	if (Ee / E0 <= 0.5 * (1 - beta) || Ee / E0 >= 0.5 * (1 + beta)) {
		std::cerr << "ERROR, Energy outside limits for PP [Lee96]! "
				<< std::endl;
		return 0.;
	} else {
		double q = E0 - Ee;
		double k = (3.0 / 4.0) * (SigmaThompson * ElectronMass * ElectronMass)
				/ (s * E0);
		double A = Ee / q + q / Ee;
		double B = E0 * (1 - beta * beta) * (1. / Ee + 1. / q);
		double C = -((1 - beta * beta) * (1 - beta * beta) * E0 * E0 / 4.0)
				* pow_integer<2>(1. / Ee + 1. / q);

		double dsigmadE = k * (A + B + C);

		return dsigmadE;
	}
}

///	 Differential cross-section for pair production for x = E/E0
double dSigmadE_PPx(double x, double beta) {

	if ((x - 0.5 * (1 - beta)) < -1 * DBL_EPSILON || x - 0.5 * (1 + beta) > DBL_EPSILON ) {
		std::cerr << "ERROR, Energy outside limits for PP [Lee96]! " << std::endl;
		std::cerr << " x = " << x << " 0.5* (1-beta) = " << 0.5 * (1 - beta) << " 0.5* (1+beta) = " << 0.5 * (1 + beta) << " beta = " << beta << std::endl;
		return 0.;
	} else {

		const double A = (x / (1. - x) + (1. - x) / x );
		const double B =  (1. / x + 1. / (1. - x) );

		return A + (1. - beta*beta) * B - (1. - beta*beta) * (1. - beta*beta) / 4. * B*B;
	}
}


/// Hold an data array to interpolate the energy distribution on 
class PPSecondariesEnergyDistribution
{
	private:
		double *_data;
		size_t _Ns;
		size_t _Nrer;
		double _s_min;
		double _s_max;
		double _dls;

	public:
		PPSecondariesEnergyDistribution(double s_min = 4. * ElectronMass * ElectronMass, double s_max =1e21,
				size_t Ns = 1000, size_t Nrer = 1000 )
		{
			if (s_min < 4.*ElectronMass*ElectronMass)
			{
				std::cerr << "Warning: Minimum COM Energy in PP Interpolation s = " << s_min << " <  (2*m_e)**2 selected. Setting to s_min = (2*m_e)**2.\n" ;
				s_min = 4.*ElectronMass*ElectronMass;
			}
			_Ns = Ns;
			_Nrer = Nrer;
			_s_min =s_min;
			_s_max = s_max;
			_data = new double[Ns*Nrer];

			_dls = (log(s_max) - log(s_min)) / (Ns);

			for (size_t i = 0; i < Ns; i++)
			{
				const double s = s_min * exp(i*_dls);
				double beta = sqrt(1. - 4. * ElectronMass*ElectronMass /s);
				
				double x0 = log((1.-beta) / 2.);
				double dx = ( log((1. + beta)/2) -  log((1.-beta) / 2.)) / (Nrer); 
				_data[i * Nrer] = exp(x0) ; 
				for (size_t j = 1; j < Nrer; j++)
				{
					double x = exp(x0 + j*dx); 
					_data[i * Nrer + j] =	dSigmadE_PPx(x, beta) + _data[i * Nrer + j - 1];
				}
			}
		}

		// returns pointer to the the integrated distribution for a given s
		double* getDistribution(double s)
		{
			size_t idx = (log(s / _s_min)) / _dls;
			double *s0 = &_data[idx * _Nrer];
			return s0;
		}

		//samples the integrated distribution and returns Eer(Ee, s)
		double sample(double E0, double eps, double theta)
		{
			double s = 2. * E0 * eps * (1-cos(theta));
			
			double *s0 = getDistribution(s); 
			double rnd = Uniform(0, 1.0) *s0[_Nrer-1];

			for (size_t i=0; i < _Nrer; i++)
			{
				if (rnd < s0[i])
				{
					double beta = sqrt(1. - 4.* ElectronMass * ElectronMass / s);

					double x0 = log((1.-beta) / 2.);
					double dx = ( log((1. + beta)/2) -  log((1.-beta) / 2.)) / (_Nrer); 
					if (Uniform(0, 1.0) < 0.5)
						return exp(x0 + (i)*dx) * E0;
					else
						return E0 * (1-exp(x0 + (i)*dx) );
				}
			}
			std::cerr << "PPSecondariesEnergyDistribution out of bounds!" << std::endl;
			std::cerr << "  s0[0] = " << s0[0] << "  s0[_Nrer-1] = " << s0[_Nrer-1] << "  rnd = " << rnd << std::endl;
			throw std::runtime_error("Grave logic error in PPSecondariesEnergyDistribution!");
		}
};



// Helper function for actual Monte Carlo sampling to avoid code-duplication
double __extractPPSecondariesEnergy(double E0, double eps, double beta)
{
	double theta = M_PI;

	static PPSecondariesEnergyDistribution interpolation;
	return interpolation.sample(E0, eps, theta);
}


double ExtractPPSecondariesEnergy(Particle &pi, Particle &pt) {
	/*!
	 Input: incident gamma Energy E0, background photon energy eps,
	 incidence angle theta.
	 Returns the energy of the produced e+ (e-)
	 */
	double E0 = pi.GetEnergy();
	double eps = pt.GetEnergy();
	double beta = pi.GetBeta();

	return __extractPPSecondariesEnergy(E0, eps, beta);
}



double ExtractPPSecondariesEnergy(Process &proc) {
	/*!
	 Input: incident gamma Energy E0, background photon energy eps,
	 incidence angle theta.
	 Returns the energy of the produced e+ (e-)
	 */

	double E0 = proc.GetIncidentParticle().GetEnergy();
	double s = proc.GetCMEnergy();
	double eps = proc.GetTargetParticle().GetEnergy();
	double theta = M_PI;
  s = ElectronMass * ElectronMass + 2 * eps * E0 * (1 - cos(theta));
	double beta = sqrt(1 - 4 * ElectronMass * ElectronMass / s);
	//	double s2 = ElectronMass * ElectronMass

	return __extractPPSecondariesEnergy(E0, eps, beta);
}


/// Hold an data array to interpolate the energy distribution on 
class ICSSecondariesEnergyDistribution
{
	private:
		double *_data;
		size_t _Ns;
		size_t _Nrer;
		double _s_min;
		double _s_max;
		double _dls;

	public:
		ICSSecondariesEnergyDistribution(double s_min = 1.01 * ElectronMass * ElectronMass /*2.6373E+11*/, double s_max =1e21,
				size_t Ns = 1000, size_t Nrer = 1000 )
		{
			// ToDo: this boundary is just an estimate
			const double l = 1.01;
			if (s_min < l * ElectronMass*ElectronMass)
			{
				std::cerr << "Warning: Minimum COM Energy in ICS Interpolation s = " << s_min << " < " << l << " m_e**2 selected. Setting to s_min = " << l << " m_e**2.\n" ;
				s_min = l * ElectronMass*ElectronMass;
			}
			_Ns = Ns;
			_Nrer = Nrer;
			_s_min =s_min;
			_s_max = s_max;
			_data = new double[Ns*Nrer];

			double theta = M_PI;

			_dls = (log(s_max) - log(s_min)) / (Ns);
			double dls_min = log(s_min);

			for (size_t i = 0; i < Ns; i++)
			{
				const double s = exp(dls_min + i*_dls);
				double beta = (s - ElectronMass * ElectronMass) / (s +
						ElectronMass * ElectronMass);

				double eer_0 = log((1-beta) / (1+beta));
				double deer =  - log((1-beta) / (1+beta)) / (Nrer);
				
				const double Ee = 1E21;
				_data[i * Nrer] = dSigmadE_ICS(Ee, Ee * exp(eer_0), s, theta); 
				for (size_t j = 1; j < Nrer; j++)
				{
					
					double Eer = Ee * exp(eer_0 + (j)*deer); 
					_data[i * Nrer + j] =	dSigmadE_ICS(Ee, Eer , s, theta) + _data[i * Nrer + j - 1];
				}
			}
		}

		// returns pointer to the the integrated distribution for a given s
		double* getDistribution(double s)
		{
			size_t idx = (log(s / _s_min)) / _dls;
			double *s0 = &_data[idx * _Nrer];
			return s0;
		}

		//samples the integrated distribution and returns Eer(Ee, s)
		double sample(double Ee, double s)
		{
			double *s0 = getDistribution(s); 
			double rnd = Uniform(0, 1.0) *s0[_Nrer-1];
			for (size_t i=0; i < _Nrer; i++)
			{
				if (rnd < s0[i])
				{
					double beta = (s - ElectronMass * ElectronMass) / (s +
								ElectronMass * ElectronMass);
					double eer_0 = log((1-beta) / (1+beta));
					double deer =  - log((1-beta) / (1+beta)) / (_Nrer );
					return exp(eer_0 + (i)*deer) * Ee; 
				}
			}
			throw std::runtime_error("Grave logic error in sampling ICSSecondariesEnergyDistribution!");	
		}
};


// Helper function for actual Monte Carlo sampling to avoid code-duplication
double __extractICSSecondaries(double Ee, double s, double theta)
{

	static ICSSecondariesEnergyDistribution interpolation;
	return interpolation.sample(Ee, s);

	//double beta = (s - ElectronMass * ElectronMass)
	//		/ (s + ElectronMass * ElectronMass);
	//bool failed = 1;

	//// reInitialization to zero..
	//double MC_Sampling_Hist[MC_SAMPLING][3];
	//for (int i = 0; i < MC_SAMPLING; i++) {
	//	for (int j = 0; j < 3; j++)
	//		MC_Sampling_Hist[i][j] = 0.;
	//}

	//double f = pow((double) (1 + beta) / (1 - beta), (double) 1. / MC_SAMPLING);
	//int cnt = 0;
	//double NormFactor = 0;

	//for (double Eer = f * ((1 - beta) / (1 + beta)) * Ee; Eer <= Ee; Eer *= f) {
	//	MC_Sampling_Hist[cnt][0] = Eer;
	//	MC_Sampling_Hist[cnt][1] = dSigmadE_ICS(Ee, Eer, s, theta);

	//	NormFactor += MC_Sampling_Hist[cnt][1];
	//	MC_Sampling_Hist[cnt][2] = NormFactor;

	//	if (MC_Sampling_Hist[cnt][1] > 0.) {
	//		cnt++;
	//	} else {
	//		break;
	//	}
	//}

	//NormFactor = (double) 1. / (double) NormFactor;

	//for (int i = 0; i < cnt; i++)
	//	MC_Sampling_Hist[i][2] *= NormFactor;

	//double rnd = 0;
	//double Eer = 0;

	//while (failed) {
	//	rnd = Uniform(0, 1.0);
	//	Eer = 0;
	//	for (int i = 0; i < cnt - 1; i++) {
	//		if (MC_Sampling_Hist[i][2] <= rnd <= MC_Sampling_Hist[i + 1][2]) {
	//			Eer = MC_Sampling_Hist[i][0];
	//			failed = 0;
	//			break;
	//		}
	//	}
	//}
	//return Eer;
}


double ExtractICSSecondariesEnergy(Particle &pi, Particle &pt) {
	/*!
	 Input: incident electron energy Ee, background photon energy eps,
	 incidence angle theta.
	 Returns the energy of the recoiled e+ (e-)
	 */
	if (::abs(pi.GetType()) != 11) {
		std::cerr << "something wrong in type ExtractICSEnergy " << std::endl;
		return 0.;
	}

	double Ee = pi.GetEnergy();
	double eps = pt.GetEnergy();
	double theta = M_PI;
	double s = 2 * Ee * eps * (1 - pi.GetBeta() * cos(cPI))
			+ pi.GetMass() * pi.GetMass();

	return __extractICSSecondaries(Ee, s, theta);
}


double ExtractICSSecondariesEnergy(Process &proc) {
	/*!
	 Input: incident electron energy Ee, background photon energy eps,
	 incidence angle theta.
	 Returns the energy of the recoiled e+ (e-)
	 */
	double Ee = proc.GetIncidentParticle().GetEnergy();
	double s = proc.GetCMEnergy();
	double theta = proc.GetInteractionAngle();
	return __extractICSSecondaries(Ee, s , theta);
}


double ExtractTPPSecondariesEnergy(Particle &pi, Particle &pt) {
	/* approximation based on A. Mastichiadis et al.,
	 Astroph. Journ. 300:178-189 (1986), eq. 30.
	 This approx is valid only for   alpha >=100
	 where alpha = p0*eps*costheta - E0*eps;
	 for our purposes, me << E0 --> p0~ E0 -->
	 alpha = E0*eps*(costheta - 1) >= 100;
	 */

	double E0 = pi.GetEnergy();
	double eps = pt.GetEnergy();
	double s = 2 * E0 * eps * (1 - pi.GetBeta() * cos(M_PI))
			+ pi.GetMass() * pi.GetMass();
	double Epp = 5.7e-1 * pow(eps / ElectronMass, -0.56) * pow(E0 / ElectronMass, 0.44) * ElectronMass;
	double Epp2 = E0
			* (1 - 1.768 * pow(s / ElectronMass / ElectronMass, -3.0 / 4.0))
			/ 2.0;
	//return the energy of each e+/e- in the pair.
	return Epp;
}


double ExtractTPPSecondariesEnergy(Process &proc) {
	/* approximation based on A. Mastichiadis et al.,
	 Astroph. Journ. 300:178-189 (1986), eq. 30.
	 This approx is valid only for   alpha >=100
	 where alpha = p0*eps*costheta - E0*eps;
	 for our purposes, me << E0 --> p0~ E0 -->
	 alpha = E0*eps*(costheta - 1) >= 100;
	 */

	double E0 = proc.GetIncidentParticle().GetEnergy();
	double eps = proc.GetTargetParticle().GetEnergy();
	double Epp = 5.7e-1 * pow(eps/ElectronMass, -0.56) * pow(E0/ElectronMass, 0.44) * ElectronMass;
	double s = proc.GetCMEnergy();
	double Epp2 = E0
			* (1 - 1.768 * pow(s / ElectronMass / ElectronMass, -3.0 / 4.0))
			/ 2.0;
	return Epp;
}


double ExtractDPPSecondariesEnergy(double E0) {
	/*
	 we use the same assumption of lee (i.e., all the energy goes equaly shared between only 1 couple of e+e-.
	 In DPPpaper has been shown that this approximation is valid within -1.5%
	 */
	if (E0 == 0)
		std::cout << "error in extracting DPP: can not be =0 " << std::endl;
	return (double) E0 / 2.0;
}

} // namespace eleca

#endif // ELECA_ENERGY_LOSS_H_

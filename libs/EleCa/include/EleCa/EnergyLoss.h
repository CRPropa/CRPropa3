#ifndef ELECA_ENERGY_LOSS_H_
#define ELECA_ENERGY_LOSS_H_

namespace eleca {

class Particle;
class Process;

double dTdZ(double z);
double betaRsh(double z);
double fLossAdiabatic(double E, double z);
double AdiabaticELoss(double z0, double z, double E0);
double MeanRateSynchrotronLoss(double E, double B);
double ESynchrotronAdiabaticLoss(double z, double E, double B);

void InitRK();
double EnergyLoss1D(double Energy, double z0, double zfin, double B);
double dSigmadE_ICS(double Ee, double Eer, double s, double theta);
double dSigmadE_PP(double Ee, double E0, double eps, double theta, double s);
double ExtractPPSecondariesEnergy(Particle &pi, Particle &pt);
double ExtractPPSecondariesEnergy(Process &proc);
double ExtractICSSecondariesEnergy(Particle &pi, Particle &pt);
double ExtractICSSecondariesEnergy(Process &proc);
double ExtractTPPSecondariesEnergy(Particle &pi, Particle &pt);
double ExtractTPPSecondariesEnergy(Process &proc);
double ExtractDPPSecondariesEnergy(double E0);

} // namespace eleca

#endif // ELECA_ENERGY_LOSS_H_

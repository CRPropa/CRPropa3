#ifndef ELECA_PROPAGATION_H
#define ELECA_PROPAGATION_H

#include "EleCa/Particle.h"
#include "EleCa/Process.h"

#include <vector>
#include <string>

namespace eleca {

class Process;
class Propagation {

private:

	double vPPle[1101];
	double vDPPle[1101];
	double vTPPle[1101];
	double vICSle[1101];
	double vEtab[1101];

	std::vector<double> BkgE, BkgA;
	std::string Bkg;
	double fEthr;
	double _dEtab;

	double magneticFieldStrength;
public:

	Propagation();

	void SetEthr(double eth) {fEthr = eth;};
	double GetEthr();
	~Propagation();

	void WriteOutput(std::ostream &out, Particle &p1,
			std::vector<Particle> &part) const;

	void ReadTables(const std::string &file);
	void InitBkgArray(const std::string &BackRad);

	double GetMeanThetaBFDeflection(double Bin, double Ein, int ptype,
			double Lin) const;
	double GetLambdaTab(const Process &proc, Process::Name procName) const;
	double ExtractMinDist(Process &proc, int type, double R, double R2,
			std::vector<double> &Etarget) const;
	std::vector<double> GetEtarget(Process &proc,
			const Particle &particle) const;
	void Propagate(Particle &curr_particle,
			std::vector<Particle> &ParticleAtMatrix,
			std::vector<Particle> &ParticleAtGround) const;
	double ExtractPhotonEnergyMC(double z, Process &proc) const;
	double ShootPhotonEnergyMC(double z) const;
	double ShootPhotonEnergyMC(double Emin, double z) const;
	void SetInitVar(std::vector<std::vector<double> > bk,
			std::vector<std::vector<double> > *le) const;

	void SetB(double B)
	{
		magneticFieldStrength = B;
	}
};

} // namespace eleca

#endif // ELECA_PROPAGATION_H

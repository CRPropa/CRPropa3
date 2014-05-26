#ifndef ELECA_PROCESS_H_
#define ELECA_PROCESS_H_

#include "Particle.h"

#include <string>
#include <iostream>

namespace eleca {

class Process {

public:

	std::string fname;
	double flambda;
	double fsmin;
	double fsmax;
	double fCMEnergy;
	double fInteractionAngle;
	double feps_inf;
	double feps_sup;
	Particle fPi;
	Particle fPt;

	std::string fback;
	double fbackdensity;

	Process();
	Process(const Process&);
	Process(Particle&, Particle&);
	Process(Particle&, Particle&, std::string);

	~Process();

	void SetName(std::string nm);
	const std::string &GetName() const;

	void SetInteractionAngle(double a);
	double GetInteractionAngle() const;

	void SetLambda(double le);
	double GetLambda() const;

	void SetLimits(double smin, double smax);
	void SetLimits(Particle& p1, std::string nameproc);
	void SetLimits();

	void SetMax(double smax);
	void SetMin(double smin);
	double GetMin() const;
	double GetMax() const;

	void SetCMEnergy(double s);

	void SetCMEnergy(Particle p1, Particle pb);

	void SetCMEnergy();
	double GetCMEnergy() const;

	void SetIncidentParticle(const Particle& p1);
	void SetTargetParticle(Particle& p1);
	const Particle &GetIncidentParticle() const;
	const Particle &GetTargetParticle() const;

	const std::string &GetBackground() const;
	void SetBackground(std::string BackRad);

private:

};

} // namespace

#endif

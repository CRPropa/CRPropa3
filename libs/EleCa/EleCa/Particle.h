#ifndef ELECA_PARTICLE_H_
#define ELECA_PARTICLE_H_

#include "Constants.h"

namespace eleca {

class Particle {
private:
	int ftype;
	double fE0ph;
	double fz0ph;
	double fbeta;
	double fmass;
	bool fIsGood;
	int fwi;
	double fB;

public:
	Particle();
	Particle(int _ft, double _fE, double fz);
//	Particle(Candidate candidate);
	~ Particle() {
	}

	bool IsGood() {
		double zmax_prop0 = 0.91425;
		double zmax_prop1 = -0.101717;
		double zmax_prop2 = 0.002855;
		double Eloc = log10(fE0ph);
		if (Eloc < 18)
			Eloc = 18;
		if (fz0ph > zmax_prop0 + zmax_prop1 * Eloc + zmax_prop2 * Eloc * Eloc)
			return 0;
		return 1;
	}

	int GetType() {
		return ftype;
	}

	void SetType(int _ft) {
		ftype = _ft;
	}

	int GetWeigth() {
		return fwi;
	}

	void SetWeigth(int _wi) {
		fwi = _wi;
	}

	double GetEnergy() {
		return fE0ph;
	}

	void SetEnergy(double _fE) {
		fE0ph = _fE;
		SetBetaAndMass();
	}

	double Getz() {
		return fz0ph;
	}

	void Setz(double _fz) {
		fz0ph = _fz;
	}

	double GetMass() {
		return fmass;
	}

	double GetBeta() {
		return fbeta;
	}

	bool GetStatus() {
		fIsGood = IsGood();
		return fIsGood;
	}

	void SetBetaAndMass() {
		if (ftype == 22) {
			fbeta = 1.0;
			fmass = 0;
		}
		if (abs(ftype) == 11) {
			fmass = ElectronMass;
			fbeta = (double) sqrt(1 - fmass * fmass / (fE0ph * fE0ph));

		}
	}

	void SetB(double B) {
		fB = B;
	}
	double GetB() {
		return fB;
	}

};

Particle::Particle(int _ft, double _fE, double _fz) {
	ftype = _ft;
	fE0ph = _fE;
	fz0ph = _fz;
	SetBetaAndMass();
	fIsGood = IsGood();
	fB = 0;
	fwi = 1;
}

Particle::Particle() {
	ftype = 22;
	fE0ph = 0;
	fz0ph = 0;
	fmass = 0;
	fbeta = 0;
	fIsGood = 0;
	fB = 0;
	fwi = 1;
}

} // namespace eleca

#endif // ELECA_PARTICLE_H_

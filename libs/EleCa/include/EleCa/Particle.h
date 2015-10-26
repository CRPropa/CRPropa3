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
  int fgeneration;

public:
	Particle();
	Particle(int _ft, double _fE, double _fz, int _fgeneration = 0);
	~ Particle();
	bool IsGood();
	int GetType() const;

	void SetType(int _ft);

	int GetWeigth() const;

	void SetWeigth(int _wi);

	double GetEnergy() const;

	void SetEnergy(double _fE);

	double Getz() const;

	void Setz(double _fz);

	double GetMass() const;

	double GetBeta() const;
//
//	bool GetStatus() {
//		fIsGood = IsGood();
//		return fIsGood;
//	}

	void SetBetaAndMass();

  int Generation();



};

} // namespace eleca

#endif // ELECA_PARTICLE_H_

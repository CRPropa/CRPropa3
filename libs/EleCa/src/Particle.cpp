#include "EleCa/Particle.h"
#include <cstdlib>

namespace eleca {

Particle::~Particle() {
}

bool Particle::IsGood() {
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

int Particle::GetType() const {
	return ftype;
}

void Particle::SetType(int _ft) {
	ftype = _ft;
}

int Particle::GetWeigth() const {
	return fwi;
}

void Particle::SetWeigth(int _wi) {
	fwi = _wi;
}

double Particle::GetEnergy() const {
	return fE0ph;
}

void Particle::SetEnergy(double _fE) {
	fE0ph = _fE;
	SetBetaAndMass();
}

double Particle::Getz() const {
	return fz0ph;
}

void Particle::Setz(double _fz) {
	fz0ph = _fz;
}

double Particle::GetMass() const {
	return fmass;
}

double Particle::GetBeta() const {
	return fbeta;
}

void Particle::SetBetaAndMass() {
	if (ftype == 22) {
		fbeta = 1.0;
		fmass = 0;
	}
	if (abs(ftype) == 11) {
		fmass = ElectronMass;
		fbeta = (double) sqrt(1 - fmass * fmass / (fE0ph * fE0ph));

	}
}

void Particle::SetB(double B) {
	fB = B;
}
double Particle::GetB() const {
	return fB;
}

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

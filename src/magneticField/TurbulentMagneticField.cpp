#include "crpropa/magneticField/TurbulentMagneticField.h"
#include "crpropa/Units.h"

namespace crpropa {

TurbulentMagneticField::TurbulentMagneticField(double Brms, double lMin,
		double lMax, double spectralIndex, int nModes) {
	setTurbulenceProperties(Brms, lMin, lMax, spectralIndex, nModes);
	initialize();
}

Vector3d TurbulentMagneticField::getField(const Vector3d &position) const {
	Vector3d b(0.);
	double a;
	for (int i = 0; i < nModes; i++) {
		Mode mode = modes[i];
		a = mode.k.dot(position) + mode.phase;
		b += mode.amplitude * (cos(a) * mode.e1 - sin(a) * mode.e2);
	}
	return b;
}

void TurbulentMagneticField::setTurbulenceProperties(double Brms, double lMin,
		double lMax, double spectralIndex, int nModes) {
	this->nModes = nModes;
	this->lMin = lMin;
	this->lMax = lMax;
	this->spectralIndex = spectralIndex;
	this->Brms = Brms;
}

void TurbulentMagneticField::initialize() {
	if (lMin >= lMax)
		throw std::runtime_error("crpropa::TurbulentMagneticField: lMin >= lMax");

	double lkMin = log10(2 * M_PI / lMax);
	double lkMax = log10(2 * M_PI / lMin);
	double dlk = (lkMax - lkMin) / (nModes - 1);
	double Lc = getCorrelationLength();

	Mode mode;
	Vector3f ek, e1, e2; // orthogonal base
	Vector3f n0(1, 1, 1); // arbitrary vector to construct orthogonal base

	double sumGk = 0;
	for (int i = 0; i < nModes; i++) {
		// construct an orthogonal base ek, e1, e2
		ek = random.randVector();
		if (ek.getAngleTo(n0) < 1e-3) {
			// ek parallel to (1,1,1)
			e1.setXYZ(-1., 1., 0);
			e2.setXYZ(1., 1., -2.);
		} else {
			// ek not parallel to (1,1,1)
			e1 = n0.cross(ek);
			e2 = ek.cross(e1);
		}
		double k = pow(10, lkMin + i * dlk);
		mode.k = ek * k;

		// amplitude (this seems to be wrong)
		double dk = k * dlk;
		double Gk = k * k * dk / (1 + pow(k * Lc, spectralIndex));
		sumGk += Gk;
		mode.amplitude = Brms * sqrt(Gk);

		// random orientation of b
		double alpha = random.rand(2 * M_PI);
		mode.e1 = e1 / e1.getR() * cos(alpha);
		mode.e2 = e2 / e2.getR() * sin(alpha);

		// random phase
		mode.phase = random.rand(2 * M_PI);
		modes.push_back(mode);
	}

	for (int i = 0; i < nModes; i++)
		modes[i].amplitude /= sumGk;
}

void TurbulentMagneticField::initialize(int seed) {
	random.seed(seed);
	initialize();
}

double TurbulentMagneticField::getPowerSpectralIndex() const {
	return spectralIndex;
}

double TurbulentMagneticField::getMinimumWavelength() const {
	return lMin;
}

double TurbulentMagneticField::getMaximumWavelength() const {
	return lMax;
}

double TurbulentMagneticField::getRMSFieldStrength() const {
	return Brms;
}

double TurbulentMagneticField::getCorrelationLength() const {
	double r = lMin / lMax;
	double a = -spectralIndex - 2;
	return lMax / 2 * (a - 1) / a * (1 - pow(r, a)) / (1 - pow(r, a - 1));
}

} // namespace crpropa

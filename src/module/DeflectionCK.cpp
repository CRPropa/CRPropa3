#include "mpc/module/DeflectionCK.h"

#include <limits>

namespace mpc {

/**
 * @class LorentzForce
 * @brief Time-derivative in SI-units of phase-point
 * (position, momentum) -> (velocity, force)
 * of a highly relativistic charged particle in magnetic field.
 */
class LorentzForce: public ExplicitRungeKutta<PhasePoint>::F {
public:
	ParticleState *particle;
	MagneticField *field;

	LorentzForce(ParticleState *particle, MagneticField *field) {
		this->particle = particle;
		this->field = field;
	}

	PhasePoint operator()(double t, const PhasePoint &v) {
		Vector3 velocity = v.b.unit() * c_light;
		Vector3 B = field->getField(v.a);
		Vector3 force = (double) particle->getChargeNumber() * eplus
				* velocity.cross(B);
		return PhasePoint(velocity, force);
	}
};

DeflectionCK::DeflectionCK(MagneticField *field, ControlType controlType,
		double tolerance) {
	erk.loadCashKarp();
	this->controlType = controlType;
	this->tolerance = tolerance;
	this->field = field;
}

DeflectionCK::~DeflectionCK() {
	delete field;
}

std::string DeflectionCK::getDescription() const {
	return "Cash-Karp Runge Kutta integration";
}

void DeflectionCK::process(Candidate *candidate,
		std::vector<Candidate *> &secondaries) {

	PhasePoint yIn(candidate->current.getPosition(),
			candidate->current.getMomentum());
	PhasePoint yOut, yErr, yScale;
	LorentzForce dydt(&candidate->current, field);
	double h = candidate->getNextStep() / c_light;
	double hTry, r;

	// phase-point to compare with error for step size control
	yScale = (yIn.abs() + dydt(0., yIn).abs() * h) * tolerance;

	do {
		hTry = h;
		erk.step(0, yIn, yOut, yErr, hTry, dydt);
		if (controlType == NoStepSizeControl) {
			// no step size control
			break;
		} else if (controlType == WorstOffender) {
			// maximum of ratio yErr(i) / yScale(i)
			r = 0;
			if (yScale.b.x() > std::numeric_limits<double>::min())
				r = std::max(r, fabs(yErr.b.x() / yScale.b.x()));
			if (yScale.b.y() > std::numeric_limits<double>::min())
				r = std::max(r, fabs(yErr.b.y() / yScale.b.y()));
			if (yScale.b.z() > std::numeric_limits<double>::min())
				r = std::max(r, fabs(yErr.b.z() / yScale.b.z()));
		} else if (controlType == RMS) {
			// RMS of ratio yErr(i) / yScale(i)
			r = 0;
			if (yScale.b.x() > std::numeric_limits<double>::min())
				r += pow(yErr.b.x() / (yScale.b.x()), 2);
			if (yScale.b.y() > std::numeric_limits<double>::min())
				r += pow(yErr.b.y() / (yScale.b.y()), 2);
			if (yScale.b.z() > std::numeric_limits<double>::min())
				r += pow(yErr.b.z() / (yScale.b.z()), 2);
			r = pow(r / 3., 0.5);
		}
		// for efficient integration try to keep r close to one
		h *= 0.95 * pow(r, -0.2);
		// limit step change
		h = std::max(h, 0.1 * hTry);
		h = std::min(h, 5 * hTry);
	} while (r > 1);

	candidate->current.setPosition(yOut.a);
	candidate->current.setDirection(yOut.b.unit());
	candidate->setCurrentStep(hTry * c_light);
	candidate->setNextStep(h * c_light);
}

} /* namespace mpc */

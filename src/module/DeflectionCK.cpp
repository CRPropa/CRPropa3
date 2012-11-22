#include "mpc/module/DeflectionCK.h"

#include <limits>
#include <sstream>
#include <stdexcept>

namespace mpc {

/**
 * @class LorentzForce
 * @brief Time-derivative of the phase-point (position, momentum) of a highly relativistic particle in magnetic field.
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
		Vector3d velocity = v.b.getUnitVector() * c_light;
		Vector3d B(0, 0, 0);
		try {
			B = field->getField(v.a);
		} catch (std::exception &e) {
			std::cerr << "mpc::LorentzForce: Exception in getField."
					<< std::endl;
			std::cerr << e.what() << std::endl;
		}
		Vector3d force = (double) particle->getChargeNumber() * eplus
				* velocity.cross(B);
		return PhasePoint(velocity, force);
	}
};

DeflectionCK::DeflectionCK(ref_ptr<MagneticField> field, double tolerance,
		double minStep, double maxStep) :
		field(field), tolerance(tolerance), minStep(minStep), maxStep(maxStep) {
	erk.loadCashKarp();
}

void DeflectionCK::setField(ref_ptr<MagneticField> f) {
	field = f;
}

void DeflectionCK::setTolerance(double t) {
	tolerance = t;
}

void DeflectionCK::setMinimumStep(double s) {
	minStep = s;
}

void DeflectionCK::setMaximumStep(double s) {
	maxStep = s;
}

std::string DeflectionCK::getDescription() const {
	std::stringstream s;
	s << "Propagation in magnetic fields using the Cash-Karp method.";
	s << " Tolerance: " << tolerance;
	s << ", Minimum Step: " << minStep / kpc << " kpc";
	s << ", Maximum Step: " << maxStep / kpc << " kpc";
	return s.str();
}

void DeflectionCK::process(Candidate *candidate) const {
	// save the new previous particle state
	candidate->previous = candidate->current;

	double step = candidate->getNextStep();
	step = std::max(step, minStep);
	step = std::min(step, maxStep);

	// rectlinear propagation for neutral particles
	if (candidate->current.getCharge() == 0) {
		Vector3d pos = candidate->current.getPosition();
		Vector3d dir = candidate->current.getDirection();
		candidate->current.setPosition(pos + dir * step);
		candidate->setCurrentStep(step);
		candidate->setNextStep(maxStep);
		return;
	}

	PhasePoint yIn(candidate->current.getPosition(),
			candidate->current.getMomentum());
	PhasePoint yOut, yErr, yScale;
	LorentzForce dydt(&candidate->current, field);
	double h = step / c_light;
	double hTry, r;

	// phase-point to compare with error for step size control
	yScale = (yIn.abs() + dydt(0., yIn).abs() * h) * tolerance;

	// try performing a steps until the relative error is less than the desired tolerance
	// or the minimum step size has been reached
	do {
		hTry = h;
		erk.step(0, yIn, yOut, yErr, hTry, dydt);

		// determine maximum of relative errors yErr(i) / yScale(i)
		r = 0;
		if (yScale.b.x > std::numeric_limits<double>::min())
			r = std::max(r, fabs(yErr.b.x / yScale.b.x));
		if (yScale.b.y > std::numeric_limits<double>::min())
			r = std::max(r, fabs(yErr.b.y / yScale.b.y));
		if (yScale.b.z > std::numeric_limits<double>::min())
			r = std::max(r, fabs(yErr.b.z / yScale.b.z));

		// change (next) step size to keep the relative error close to the tolerance
		h *= 0.95 * pow(r, -0.2);
		h = std::max(h, 0.1 * hTry);
		h = std::min(h, 5 * hTry);
	} while (r > 1 && h > minStep);

	candidate->current.setPosition(yOut.a);
	candidate->current.setDirection(yOut.b.getUnitVector());
	candidate->setCurrentStep(hTry * c_light);
	candidate->setNextStep(h * c_light);
}

} /* namespace mpc */

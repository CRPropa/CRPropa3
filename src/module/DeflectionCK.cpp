#include "mpc/module/DeflectionCK.h"

#include <limits>
#include <sstream>

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
		Vector3d velocity = v.b / v.b.mag() * c_light;
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

DeflectionCK::DeflectionCK(MagneticField *field, double tolerance,
		double minimumStep) {
	this->field = field;
	this->tolerance = tolerance;
	this->minimumStep = minimumStep;
	erk.loadCashKarp();
}

std::string DeflectionCK::getDescription() const {
	std::stringstream sstr;
	sstr
			<< "Propagation in magnetic fields using the Cash-Karp method. Tolerance: "
			<< tolerance << ", Minimum Step: " << minimumStep / kpc << " kpc";
	return sstr.str();
}

void DeflectionCK::process(Candidate *candidate) const {
	// rectlinear propagation in case of no charge
	if (candidate->current.getChargeNumber() == 0) {
		double nextStep = std::max(minimumStep, candidate->getNextStep());
		Vector3d pos = candidate->current.getPosition();
		Vector3d dir = candidate->current.getDirection();
		candidate->current.setPosition(pos + dir * nextStep);
		candidate->setCurrentStep(nextStep);
		candidate->setNextStep(nextStep * 5);
		return;
	}

	PhasePoint yIn(candidate->current.getPosition(),
			candidate->current.getMomentum());
	PhasePoint yOut, yErr, yScale;
	LorentzForce dydt(&candidate->current, field);
	double h = std::max(candidate->getNextStep(), minimumStep) / c_light;
	double hTry, r;

	// phase-point to compare with error for step size control
	yScale = (yIn.abs() + dydt(0., yIn).abs() * h) * tolerance;

	do {
		hTry = h;
		erk.step(0, yIn, yOut, yErr, hTry, dydt);

		// maximum of ratio yErr(i) / yScale(i)
		r = 0;
		if (yScale.b.x > std::numeric_limits<double>::min())
			r = std::max(r, fabs(yErr.b.x / yScale.b.x));
		if (yScale.b.y > std::numeric_limits<double>::min())
			r = std::max(r, fabs(yErr.b.y / yScale.b.y));
		if (yScale.b.z > std::numeric_limits<double>::min())
			r = std::max(r, fabs(yErr.b.z / yScale.b.z));

		// for efficient integration try to keep r close to one
		h *= 0.95 * pow(r, -0.2);

		// limit step change
		h = std::max(h, 0.1 * hTry);
		h = std::min(h, 5 * hTry);
	} while (r > 1 && h >= minimumStep);

	candidate->current.setPosition(yOut.a);
	candidate->current.setDirection(yOut.b / yOut.b.mag());
	candidate->setCurrentStep(hTry * c_light);
	candidate->setNextStep(h * c_light);
}

void DeflectionCK::updateSimulationVolume(const Vector3d &origin, double size) {
	field->updateSimulationVolume(origin, size);
}

} /* namespace mpc */

#ifndef DEFLECTION_H_
#define DEFLECTION_H_

#include "Particle.h"
#include "ExplicitRungeKutta.h"
#include "ThreeVector.h"
#include "MagneticField.h"
#include "PhasePoint.h"
#include <limits>

/**
 * @class LorentzF
 * @brief Time-derivative in SI-units of phase-point
 * (position, momentum) -> (velocity, force)
 * of a highly relativistic charged particle in magnetic field.
 */
class LorentzForce: public ExplicitRungeKutta<PhasePoint>::F {
public:
	Particle *particle;
	MagneticField *field;

	LorentzForce(Particle *particle, MagneticField *field) {
		this->particle = particle;
		this->field = field;
	}

	PhasePoint operator()(double t, const PhasePoint &v) {
		Hep3Vector velocity = v.b.unit() * c_light;
		Hep3Vector B = field->getField(v.a);
		Hep3Vector force = (double) particle->getCharge() * velocity.cross(B);
		return PhasePoint(velocity, force);
	}
};

/**
 * @class Deflection
 * @brief Magnetic deflection in 3D using a the Cash-Karp Runge-Kutta method
 * propagates the particle by a step particle.getNextStep() or smaller.
 * The step size control tries to keep the error close to, but smaller than the maxError
 */
class DeflectionCK {
public:
	double tolerance;
	ExplicitRungeKutta<PhasePoint> erk;
	enum ControlType {
		NoStepSizeControl, WorstOffender, RMS
	};
	ControlType controlType;

	DeflectionCK(ControlType controlType, double tolerance) {
		erk.loadCashKarp();
		this->controlType = controlType;
		this->tolerance = tolerance;
	}

	void apply(Particle &particle, MagneticField &field) {
		LorentzForce dydt(&particle, &field);
		PhasePoint yIn(particle.getPosition(), particle.getMomentum());
		PhasePoint yOut, yErr, yScale;
		double nextStep = particle.getNextStep() / c_light;
		double step, r;

		yScale = (yIn.abs() + dydt(0., yIn).abs() * nextStep) * tolerance;

		do {
			step = nextStep;
			erk.step(0, yIn, yOut, yErr, step, dydt);

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
			nextStep *= 0.95 * pow(r, -0.2);
			// limit step change
			nextStep = std::max(nextStep, 0.1 * step);
			nextStep = std::min(nextStep, 5 * step);
		} while (r > 1);

		particle.setPosition(yOut.a);
		particle.setDirection(yOut.b.unit());
		particle.setStep(step * c_light);
		particle.setNextStep(nextStep * c_light);
	}

};

#endif /* DEFLECTION_H_ */


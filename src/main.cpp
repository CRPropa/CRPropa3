#include "mpc/Candidate.h"
#include "mpc/DeflectionCK.h"
#include "mpc/BreakCondition.h"
#include "mpc/Propagator.h"
#include "mpc/GlutDisplay.h"
#include "mpc/ParticleState.h"
#include "mpc/SharedPointer.h"

using mpc::Mpc;

class SimpleField: public mpc::Feature {
	mpc::HomogeneousMagneticField field;
	mpc::DeflectionCK deflection;
public:
	SimpleField() :
			field(mpc::Vector3(0., 0., 1e-13)), deflection(
					mpc::DeflectionCK::WorstOffender, 5e-5) {

	}

	void apply(mpc::Candidate &candidate, size_t priority) {
		deflection.apply(candidate, field);
	}

	std::string description() const {
		return "SimpleField";
	}
};

int main() {
	mpc::Propagator propagator;

	propagator.add(mpc::Priority::AfterIntegration,
			new mpc::MaximumTrajectoryLength(100 * Mpc));
	propagator.add(mpc::Priority::Integration, new SimpleField);
	propagator.add(mpc::Priority::AfterCommit, new mpc::GlutDisplay(5.));

	propagator.print();

	mpc::ParticleState initial;
	initial.setPosition(mpc::Vector3(-1.08, 0., 0.) * Mpc);
	initial.setChargeNumber(1);
	initial.setMass(1 * mpc::amu);
	initial.setDirection(mpc::Vector3(0., 1., 0.));
	initial.setEnergy(1 * mpc::EeV);

	mpc::Candidate candidate;
	candidate.next = initial;
	candidate.initial = initial;
	candidate.setNextStepSize(0.05 * mpc::Mpc);

	propagator.apply(candidate);

	propagator.clear();

	propagator.add(mpc::Priority::Integration, new SimpleField);
	candidate.next = initial;
	candidate.initial = initial;
	candidate.setNextStepSize(0.05 * mpc::Mpc);
	propagator.print();
	propagator.apply(candidate);

	return 0;
}


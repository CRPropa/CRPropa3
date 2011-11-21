#include "mpc/Candidate.h"
#include "mpc/DeflectionCK.h"
#include "mpc/BreakCondition.h"
#include "mpc/Propagator.h"
#include "mpc/BreakCondition.h"
#include "mpc/GlutDisplay.h"

using mpc::Mpc;

class SimpleField: public mpc::Feature {
	mpc::HomogeneousMagneticField field;
	mpc::DeflectionCK deflection;
public:
	SimpleField() :
			field(mpc::Vector3(0., 0., 1e-13)), deflection(
					mpc::DeflectionCK::WorstOffender, 5e-5) {

	}

	virtual void apply(mpc::Candidate &candidate, size_t priority) {
		deflection.apply(candidate, field);
	}
};

int main() {
	mpc::Propagator propa;

	mpc::MaximumTrajectoryLength maxLen(100 * Mpc);
	propa.add(mpc::Priority::AfterIntegration, &maxLen);

	SimpleField f;
	propa.add(mpc::Priority::Integration, &f);

	mpc::GlutDisplay display(5.);
	propa.add(mpc::Priority::AfterCommit, &display);

	propa.print();

	mpc::ParticleState initial;
	initial.setPosition(mpc::Vector3(-1.08, 0., 0.) * Mpc);
	initial.setChargeNumber(1);
	initial.setMass(1 * mpc::amu);
	initial.setDirection(mpc::Vector3(0., 1., 0.));
	initial.setEnergy(1 * mpc::EeV);

	mpc::Candidate candidate;
	candidate.next = initial;
	candidate.initial = initial;
	candidate.setNextStep(0.05 * mpc::Mpc);

	propa.apply(candidate);

#if 0

	mpc::MaximumTrajectoryLength maxTrajLength(10 * Mpc);

	for (int i = 0; i < 800; i++) {
//		std::cout << "step:  " << particle.getStep()
//			<< ", position (Mpc):  " << particle.getPositionMpc()
//			<< ", momentum (EeV/c):  " << particle.getMomentum() / EeV
//			<< ", check  " << particle.getDirection().mag() << std::endl;
		std::cout << candidate.getNextStep() / Mpc << ", ";
		std::cout << candidate.next.getPosition().x() / Mpc << ","
		<< candidate.next.getPosition().y() / Mpc << ","
		<< candidate.next.getPosition().z() / Mpc << std::endl;
		deflection.apply(candidate, field);
		maxTrajLength.apply(candidate);
	}
#endif

	return 0;
}


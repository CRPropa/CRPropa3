#include "mpc/Particle.h"
#include "mpc/DeflectionCK.h"
#include "mpc/BreakCondition.h"
#include "mpc/ThreeVector.h"

class Feature {

};

class Propagation {
public:
};

using mpc::Mpc;

int main() {
	Propagation propa;

	mpc::Particle initial;
	initial.setPosition(mpc::Hep3Vector(-1.08, 0., 0.) * Mpc);
	initial.setChargeNumber(1);
	initial.setMass(1 * mpc::amu);
	initial.setDirection(mpc::Hep3Vector(0., 1., 0.));
	initial.setEnergy(1 * mpc::EeV);

	mpc::Candidate candidate;
	candidate.current = initial;
	candidate.initial = initial;
	candidate.setNextStep(0.05 * mpc::Mpc);

	mpc::HomogeneousMagneticField field(mpc::Hep3Vector(0., 0., 1e-13)); // 1nG perpendicular to p
	mpc::DeflectionCK deflection(mpc::DeflectionCK::WorstOffender, 5e-5);

	mpc::MaximumTrajectoryLength maxTrajLength(10 * Mpc);

	for (int i = 0; i < 800; i++) {
//		std::cout << "step:  " << particle.getStep()
//			<< ", position (Mpc):  " << particle.getPositionMpc()
//			<< ", momentum (EeV/c):  " << particle.getMomentum() / EeV
//			<< ", check  " << particle.getDirection().mag() << std::endl;
		std::cout << candidate.getNextStep() / Mpc << ", ";
		std::cout << candidate.current.getPosition().x() / Mpc << ","
				<< candidate.current.getPosition().y() / Mpc << ","
				<< candidate.current.getPosition().z() / Mpc << std::endl;
		deflection.apply(candidate, field);
		maxTrajLength.apply(candidate);
	}

	return 0;
}


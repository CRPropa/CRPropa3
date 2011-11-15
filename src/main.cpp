#include "mpc/Particle.h"
#include "mpc/DeflectionCK.h"
#include "mpc/BreakCondition.h"
#include "mpc/ThreeVector.h"

class Feature {

};

class Propagation {
public:
};

int main() {
	Propagation propa;

	mpc::Particle initial;
	initial.setPositionMpc(mpc::Hep3Vector(-1.08, 0., 0.));
	initial.setChargeNumber(1);
	initial.setMass(1 * mpc::amu);
	initial.setDirection(mpc::Hep3Vector(0., 1., 0.));
	initial.setEnergyEeV(1);

	mpc::Candidate candidate;
	candidate.current = initial;
	candidate.initial = initial;
	candidate.setNextStepMpc(0.05);

	mpc::HomogeneousMagneticField field(mpc::Hep3Vector(0., 0., 1e-13)); // 1nG perpendicular to p
	mpc::DeflectionCK deflection(mpc::DeflectionCK::WorstOffender, 5e-5);

	mpc::MaximumTrajectoryLength maxTrajLength(10 * mpc::Mpc);

	for (int i = 0; i < 800; i++) {
//		std::cout << "step:  " << particle.getStep()
//			<< ", position (Mpc):  " << particle.getPositionMpc()
//			<< ", momentum (EeV/c):  " << particle.getMomentum() / EeV
//			<< ", check  " << particle.getDirection().mag() << std::endl;
		std::cout << candidate.getNextStepMpc() << ", ";
		std::cout << candidate.current.getPositionMpc().x() << ","
				<< candidate.current.getPositionMpc().y() << ","
				<< candidate.current.getPositionMpc().z() << std::endl;
		deflection.apply(candidate, field);
		maxTrajLength.apply(candidate);
	}

	return 0;
}


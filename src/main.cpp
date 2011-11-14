
#include "Particle.h"
#include "DeflectionCK.h"
#include "BreakCondition.h"

int main() {

	Particle particle;
	particle.setChargeMassNumber(1,1);
	particle.setStepMpc(0.05);
	particle.setNextStepMpc(0.05);
	particle.setPositionMpc(Hep3Vector(-1.08,0.,0.));
	particle.setDirection(Hep3Vector(0.,1.,0.));
	particle.setEnergyEeV(1);

	HomogeneousMagneticField field(Hep3Vector(0., 0., 1e-13)); // 1nG perpendicular to p
	DeflectionCK deflection(DeflectionCK::WorstOffender, 5e-5);

	MaximumTrajectoryLength maxTrajLength(10*Mpc);

	for (int i = 0; i < 800; i++) {
//		std::cout << "step:  " << particle.getStep()
//			<< ", position (Mpc):  " << particle.getPositionMpc()
//			<< ", momentum (EeV/c):  " << particle.getMomentum() / EeV
//			<< ", check  " << particle.getDirection().mag() << std::endl;
		std::cout << particle.getStepMpc() << ", ";
		std::cout << particle.getPositionMpc().x() << ","
				  << particle.getPositionMpc().y() << ","
				  << particle.getPositionMpc().z() << std::endl;
		deflection.apply(particle, field);
		maxTrajLength.apply(particle);
	}

	return 0;
}



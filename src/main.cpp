#include "mpc/ModuleList.h"
#include "mpc/module/SimplePropagation.h"
#include "mpc/module/DeflectionCK.h"
#include "mpc/module/BreakCondition.h"
#include "mpc/module/Output.h"
#include "mpc/module/ElectronPairProduction.h"
#include "mpc/module/PhotoPionProduction.h"
#include "mpc/module/PhotoDisintegration.h"
#include "mpc/module/NuclearDecay.h"
#include "mpc/magneticField/UniformMagneticField.h"
#include "mpc/magneticField/TurbulentMagneticField.h"

using namespace mpc;

int main(int argc, char **argv) {
	ModuleList modules;

	// propagation --------------------------------------------------------
	TurbulentMagneticField *bField = new TurbulentMagneticField(Vector3d(0, 0, 0), 64, 1., 2., 8., 1e-12, -11. / 3.);
	bField->initialize();
	modules.add(new DeflectionCK(bField));

	// interactions -------------------------------------------------------
	modules.add(new NuclearDecay());
	modules.add(new PhotoDisintegration());
	modules.add(new ElectronPairProduction());
	modules.add(new PhotoPionProduction());

	// break conditions ---------------------------------------------------
	modules.add(new MaximumTrajectoryLength(50 * Mpc));

	// output -------------------------------------------------------------
	modules.add(new ShellOutput());

	std::cout << modules << std::endl;

	ParticleState initial;
	initial.setId(getNucleusId(56, 26));
	initial.setEnergy(100 * EeV);
	initial.setPosition(Vector3d(0, 0, 0));
	initial.setDirection(Vector3d(1, 0, 0));

	ref_ptr<Candidate> candidate = new Candidate(initial);
	modules.process(candidate);

	std::cout << "done" << std::endl;

	return 0;
}

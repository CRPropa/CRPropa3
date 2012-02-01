#include "mpc/XMLImport.h"
#include "mpc/ModuleChain.h"
#include "mpc/Vector3.h"
#include "mpc/Units.h"
#include "mpc/module/DeflectionCK.h"
#include "mpc/module/BreakCondition.h"
#include "mpc/module/GlutDisplay.h"
#include "mpc/module/Output.h"
#include "mpc/module/ElectronPairProduction.h"
#include "mpc/module/PhotoPionProduction.h"
#include "mpc/module/PhotoDisintegration.h"
#include "mpc/module/NuclearDecay.h"
#include "mpc/magneticfield/TurbulentMagneticField.h"

using namespace mpc;

int main(int argc, char **argv) {
	ModuleChain chain;

	if (argc > 1) {
		XMLImport import(&chain);
		import.import(argv[1]);
	} else {
		// propagation --------------------------------------------------------
		HomogeneousMagneticField *field = new HomogeneousMagneticField(
				Vector3(0., 0., 1e-20));
//		TurbulentMagneticField field(Vector3(0, 0, 0) * Mpc, 64, 100 * kpc,
//				1. * nG, -11. / 3., 200 * kpc, 800 * kpc);
//		field.initialize();
		chain.add(new DeflectionCK(field, DeflectionCK::WorstOffender, 5e-5),
				25);

		// interactions -------------------------------------------------------
//		chain.add(new NuclearDecay(), 30);
//		chain.add(new PhotoDisintegration(), 31);
//		chain.add(new ElectronPairProduction(ElectronPairProduction::CMB), 32);
//		chain.add(new PhotoPionProduction(PhotoPionProduction::CMBIR), 33);

		// break conditions ---------------------------------------------------
//		chain.add(new MinimumEnergy(5 * EeV), 50);
//		chain.add(new MaximumTrajectoryLength(100 * Mpc), 51);
//		chain.add(new LargeObserverSphere(9 * Mpc, Vector3(0, 0, 0) * Mpc), 52);
		chain.add(new SmallObserverSphere(1 * Mpc, Vector3(5, 5, 0) * Mpc), 53);

		// output -------------------------------------------------------------
		chain.add(new ShellOutput(), 79);
//		chain.add(new TrajectoryOutput("trajectories.csv"), 80);
//		chain.add(new FinishedOutput("finished.txt"), 80);
		chain.add(new GlutDisplay(), 80);
	}

	std::cout << chain << std::endl;

	ParticleState initial;
	initial.setId(getNucleusId(56, 26));
	initial.setEnergy(100 * EeV);
	initial.setPosition(Vector3(0., 1., 0.) * Mpc);
	initial.setDirection(Vector3(1., 1., 0.));

	Candidate candidate;
	candidate.current = initial;
	candidate.initial = initial;
	candidate.setNextStep(0.01 * Mpc);

	std::vector<Candidate *> candidates;
	candidates.push_back(&candidate);

	chain.process(candidates);

	return 0;
}

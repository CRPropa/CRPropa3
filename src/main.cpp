#include "mpc/XMLImport.h"
#include "mpc/ModuleChain.h"
#include "mpc/Candidate.h"
#include "mpc/ParticleState.h"
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
		// propagation
		HomogeneousMagneticField *field = new HomogeneousMagneticField(
				Vector3(0., 0., 1e-12));
//		TurbulentMagneticField field(Vector3(0, 0, 0) * Mpc, 64, 100 * kpc,
//				1. * nG, -11. / 3., 200 * kpc, 800 * kpc);
//		field.initialize();
		chain.add(new DeflectionCK(field, DeflectionCK::WorstOffender, 5e-5),
				25);

		// interactions
//		chain.add(new NuclearDecay(), 30);
//		chain.add(new PhotoDisintegration(), 31);
//		chain.add(new ElectronPairProduction(ElectronPairProduction::CMB), 32);
		chain.add(new PhotoPionProduction(PhotoPionProduction::CMBIR), 33);

		// break conditions
		chain.add(new MinimumEnergy(5 * EeV), 50);

		// output
		chain.add(new ShellOutput(), 79);
//		chain.add(new TrajectoryOutput("trajectories.csv"), 80);
//		chain.add(new FinishedOutput("finished.txt"), 80);
		chain.add(new GlutDisplay(), 80);
	}

	std::cout << chain << std::endl;

	ParticleState initial;
	initial.setId(1001001000);
	initial.setEnergy(50 * EeV);
	initial.setPosition(Vector3(-1.08, 0., 0.) * Mpc);
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

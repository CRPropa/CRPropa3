#include "mpc/XMLImport.h"
#include "mpc/ModuleChain.h"
#include "mpc/Vector3.h"
#include "mpc/Units.h"
#include "mpc/module/SimplePropagation.h"
#include "mpc/module/DeflectionCK.h"
#include "mpc/module/BreakCondition.h"
#include "mpc/module/GlutDisplay.h"
#include "mpc/module/Output.h"
#include "mpc/module/ElectronPairProduction.h"
#include "mpc/module/PhotoPionProduction.h"
#include "mpc/module/PhotoDisintegration.h"
#include "mpc/module/NuclearDecay.h"
#include "mpc/magneticField/uniformMagneticField.h"
#include "mpc/magneticField/turbulentMagneticFieldGrid.h"
#include "mpc/magneticField/sphMagneticField.h"

using namespace mpc;

int main(int argc, char **argv) {
	ModuleChain chain;

//	if (argc > 1) {
//		XMLImport import(&chain);
//		import.import(argv[1]);
//	} else {

// propagation --------------------------------------------------------
//	UniformMagneticField *field = new UniformMagneticField(Vector3(0., 0., 1e-20));
//	SPHMagneticField *field = new SPHMagneticField(Vector3(119717, 221166, 133061) * kpc, 3 * Mpc, 50);
//	field->gadgetField->load("test/coma-0.7.raw");
//	TurbulentMagneticFieldGrid *field = new TurbulentMagneticFieldGrid(Vector3(0, 0, 0), 64, 0.1 * Mpc, 1e-12, 1, 8, -11. / 3., 10);
//	chain.add(1, new DeflectionCK(field, DeflectionCK::WorstOffender, 1e-4));
	chain.add(1, new SimplePropagation);

	// interactions -------------------------------------------------------
	chain.add(10, new NuclearDecay());
	chain.add(11, new PhotoDisintegration());
	chain.add(12, new ElectronPairProduction(ElectronPairProduction::CMB));
	chain.add(13, new PhotoPionProduction(PhotoPionProduction::CMBIR));

	// break conditions ---------------------------------------------------
	chain.add(20, new MaximumTrajectoryLength(50 * Mpc));

	// output -------------------------------------------------------------
//	chain.add(79, new ShellOutput());
//	chain.add(80, new GlutDisplay());
	TrajectoryOutput *output = new TrajectoryOutput("26-56-50-10.txt");
	chain.add(99, output);

	std::cout << chain << std::endl;

	ParticleState initial;
	initial.setId(getNucleusId(56, 26));
	initial.setEnergy(100 * EeV);
	initial.setPosition(Vector3(0, 0, 0));
	initial.setDirection(Vector3(1, 0, 0));

	std::vector<ref_ptr<Candidate> > candidates;
	for (size_t i = 0; i < 1000; i++) {
		Candidate *candidate = new Candidate;
		candidate->current = initial;
		candidate->initial = initial;
		candidates.push_back(candidate);
	}

	chain.process(candidates, true);
	delete output;

	return 0;
}

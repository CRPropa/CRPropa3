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
//	TurbulentMagneticFieldGrid *field = new TurbulentMagneticFieldGrid(Vector3(0, 0, 0), 64, 1., 2., 8., 1e-12, -11. / 3.);
//	field.initialize();
//	chain.add(1, new DeflectionCK(&field, DeflectionCK::WorstOffender, 1e-4));
	chain.add(1, new SimplePropagation);

	// interactions -------------------------------------------------------
	chain.add(10, new NuclearDecay());
	chain.add(11, new PhotoDisintegration());
	chain.add(12, new ElectronPairProduction(ElectronPairProduction::CMB));
	chain.add(13, new PhotoPionProduction(PhotoPionProduction::CMBIR));

	// break conditions ---------------------------------------------------
	chain.add(20, new MaximumTrajectoryLength(50 * Mpc));

	// output -------------------------------------------------------------
	chain.add(79, new ShellOutput());
	//chain.add(80, new GlutDisplay());

	std::cout << chain << std::endl;

	ParticleState initial;
	initial.setId(getNucleusId(56, 26));
	initial.setEnergy(100 * EeV);
	initial.setPosition(Vector3(0, 0, 0));
	initial.setDirection(Vector3(1, 0, 0));

	ref_ptr<Candidate> candidate = new Candidate(initial);

	chain.process(candidate);

	std::cout << "done" << std::endl;

	return 0;
}

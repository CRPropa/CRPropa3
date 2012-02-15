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
//	chain.add(new DeflectionCK(field, DeflectionCK::WorstOffender, 1e-4), 25);
	chain.add(new SimplePropagation(10), 25);

	// interactions -------------------------------------------------------
	chain.add(new NuclearDecay(), 30);
	chain.add(new PhotoDisintegration(), 31);
	chain.add(new ElectronPairProduction(ElectronPairProduction::CMB), 32);
	chain.add(new PhotoPionProduction(PhotoPionProduction::CMBIR), 33);

	// break conditions ---------------------------------------------------
//	chain.add(new MinimumEnergy(5 * EeV), 50);
//	chain.add(new MaximumTrajectoryLength(50 * Mpc), 51);
	chain.add(new SphericalBoundary(Vector3(0, 0, 0) * Mpc, 9 * Mpc, 1 * Mpc, Candidate::Detected), 52);
//	chain.add(new SmallObserverSphere(Vector3(0, 0, 0) * Mpc, 1 * Mpc), 53);

	// output -------------------------------------------------------------
	chain.add(new ShellOutput(), 79);
//	chain.add(new TrajectoryOutput("trajectories.csv"), 80);
//	chain.add(new GlutDisplay(), 80);
	ChargeMassEngergyOutput *output = new ChargeMassEngergyOutput("26-56-50-10.txt", Candidate::Detected);
	chain.add(output, 100);

	std::cout << chain << std::endl;

	ParticleState initial;
	initial.setId(getNucleusId(56, 26));
	initial.setEnergy(50 * EeV);
	initial.setPosition(Vector3(0, 0, 0));
	initial.setDirection(Vector3(1, 0, 0));

	std::vector<Candidate *> candidates;
	for (size_t i=0; i<1000; i++) {
		Candidate *candidate = new Candidate;
		candidate->current = initial;
		candidate->initial = initial;
		candidate->setNextStep(0.01 * Mpc);
		candidates.push_back(candidate);
	}

	chain.process(candidates);
	delete output;

	return 0;
}

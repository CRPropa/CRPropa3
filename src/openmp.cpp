#include "mpc/magneticField/turbulentMagneticFieldGrid.h"
#include "mpc/magneticField/uniformMagneticField.h"
#include "mpc/module/DeflectionCK.h"
#include "mpc/module/NuclearDecay.h"
#include "mpc/module/PhotoDisintegration.h"
#include "mpc/module/ElectronPairProduction.h"
#include "mpc/module/PhotoPionProduction.h"
#include "mpc/module/BreakCondition.h"
#include "mpc/module/Output.h"
#include "mpc/ModuleList.h"

#include <omp.h>

using namespace mpc;
using namespace std;

int main(int argc, char **argv) {

//	ref_ptr<TurbulentMagneticFieldGrid> field = new TurbulentMagneticFieldGrid(
//			Vector3(0, 0, 0), 64, 1., 2., 8., 1e-12, -11. / 3.);
	ref_ptr<UniformMagneticField> field = new UniformMagneticField(
			Vector3(0, 0, 1e-12));

	ModuleList modules;
	modules.add(new DeflectionCK(field));
	modules.add(new NuclearDecay());
	modules.add(new PhotoDisintegration());
	modules.add(new ElectronPairProduction(CMBIR));
	modules.add(new PhotoPionProduction(CMBIR));
	modules.add(new MaximumTrajectoryLength(50 * Mpc));
	modules.add(new MinimumEnergy(5 * EeV));

	cout << modules << endl;

	ParticleState initial;
	initial.setId(getNucleusId(56, 26));
	initial.setEnergy(100 * EeV);
	initial.setPosition(Vector3(0, 0, 0));
	initial.setDirection(Vector3(1, 0, 0));

#pragma omp parallel for
	for (size_t i = 0; i < 1000; i++) {

		ref_ptr<Candidate> candidate = new Candidate(initial);
		modules.run(candidate, true);
	}

	cout << "done" << endl;

	return 0;
}

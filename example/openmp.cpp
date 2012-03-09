#include "mpc/magneticField/turbulentMagneticFieldGrid.h"
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

	ref_ptr<TurbulentMagneticFieldGrid> field = new TurbulentMagneticFieldGrid(
			Vector3(0, 0, 0), 64, 1., 2., 8., 1e-12, -11. / 3.);

	ModuleList modules;
	modules.push_back(new DeflectionCK(field));
	modules.push_back(new NuclearDecay());
	modules.push_back(new PhotoDisintegration());
	modules.push_back(new ElectronPairProduction(ElectronPairProduction::CMB));
	modules.push_back(new PhotoPionProduction(PhotoPionProduction::CMBIR));
	modules.push_back(new MaximumTrajectoryLength(50 * Mpc));
	modules.push_back(new MinimumEnergy(5 * EeV));
	//modules.push_back(new ShellOutput);

	cout << modules << endl;

#pragma omp parallel  for
	for (size_t i = 0; i < 10000; i++) {
#pragma omp critical
		cout << "[" << omp_get_thread_num() << "] " << i << endl;
		ParticleState initial;
		initial.setId(getNucleusId(56, 26));
		initial.setEnergy(100 * EeV);
		initial.setPosition(Vector3(0, 0, 0));
		initial.setDirection(Vector3(1, 0, 0));

		ref_ptr<Candidate> candidate = new Candidate(initial);

		modules.run(candidate, true);
	}

	cout << "done" << endl;

	return 0;
}

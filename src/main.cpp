#include "mpc/Candidate.h"
#include "mpc/DeflectionCK.h"
#include "mpc/BreakCondition.h"
#include "mpc/ModuleChain.h"
#include "mpc/GlutDisplay.h"
#include "mpc/ParticleState.h"
#include "mpc/magneticfield/TurbulentMagneticField.h"
#include "mpc/magneticfield/MagneticFieldRing.h"
#include "mpc/Units.h"
#include "mpc/Output.h"
#include "mpc/ModuleChainImport.h"
#include "mpc/interaction/ElectronPairProduction.h"
#include "mpc/interaction/Decay.h"

#include <iostream>

using namespace mpc;

int main() {
	ModuleChain chain;

	HomogeneousMagneticField field(Vector3(0., 0., 1e-12));
//	TurbulentMagneticField field(Vector3(0, 0, 0) * Mpc, 64, 100 * kpc, 1. * nG,
//			-11. / 3., 200 * kpc, 800 * kpc);
//	field.initialize();

	chain.add(new DeflectionCK(&field, DeflectionCK::WorstOffender, 5e-5), 25);
	chain.add(new Decay(), 30);
//	chain.add(new ElectronPairProduction, 31);
	chain.add(new GlutDisplay(), 80);

	std::cout << chain << std::endl;

	ParticleState initial;
	initial.setId(1005002000);
	initial.setEnergy(100 * EeV);
	initial.setPosition(Vector3(-1.08, 0., 0.) * Mpc);
	initial.setDirection(Vector3(1., 1., 0.));

	Candidate candidate;
	candidate.current = initial;
	candidate.initial = initial;
	candidate.setNextStep(0.01 * Mpc);

	std::vector<Candidate *> event;
	event.push_back(&candidate);

	chain.process(event);

	return 0;
}

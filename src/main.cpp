#include "mpc/Candidate.h"
#include "mpc/DeflectionCK.h"
#include "mpc/BreakCondition.h"
#include "mpc/ModuleChain.h"
#include "mpc/GlutDisplay.h"
#include "mpc/ParticleState.h"
#include "mpc/TurbulentMagneticField.h"
#include "mpc/Units.h"
#include "mpc/Output.h"
#include "mpc/ModuleChainImport.h"

#include <iostream>

using namespace mpc;

int main() {
	ModuleChain chain;

	// magnetic field
	TurbulentMagneticField field(Vector3(0, 0, 0) * Mpc, 64, 100 * kpc, 1. * nG,
			-11. / 3., 2., 8.);
	field.initialize();

	// deflection
	DeflectionCK deflection(&field, DeflectionCK::WorstOffender, 5e-5);
	chain.add(Priority::Integration, &deflection);

	// maximum trajectory length
	MaximumTrajectoryLength trajLength(100 * Mpc);
	chain.add(mpc::Priority::AfterIntegration, &trajLength);

//	chain.add(Priority::AfterCommit, new GlutDisplay());
//	chain.add(Priority::AfterCommit, new CandidateOutput());

	import(chain, "example.xml");

	std::cout << chain << std::endl;

	ParticleState initial;
//	initial.setPosition(Vector3(-1.08, 0., 0.) * Mpc);
	initial.setPosition(Vector3(2, 2, 2) * Mpc);
	initial.setChargeNumber(1);
	initial.setMass(1 * amu);
	initial.setDirection(Vector3(1., 1., 1.));
	initial.setEnergy(10 * EeV);

	Candidate candidate;
	candidate.current = initial;
	candidate.initial = initial;
	candidate.setNextStep(0.05 * Mpc);

	std::vector<Candidate *> event;
	event.push_back(&candidate);

	chain.process(event);

	return 0;
}

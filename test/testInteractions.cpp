#include "gtest/gtest.h"
#include "mpc/Candidate.h"
#include "mpc/interaction/ElectronPairProduction.h"
#include <iostream>

namespace mpc {

TEST(testElectronPairProduction, EnergyDecreasing) {
	// test if energy loss occurs for protons with energies from 1e15 - 1e23 eV
	Candidate candidate;
	candidate.setCurrentStep(1 * Mpc);
	candidate.current.setId(1000010010); // proton
	std::vector<Candidate *> secondaries;

	ElectronPairProduction epp1(ElectronPairProduction::CMB);
	for (int i = 0; i<80; i++) {
		double E = pow(10, 15+i*0.1) * eV;
		candidate.current.setEnergy(E);
		epp1.process(&candidate, secondaries);
		EXPECT_TRUE(candidate.current.getEnergy() <= E);
	}

	ElectronPairProduction epp2(ElectronPairProduction::IR);
	for (int i = 0; i<80; i++) {
		double E = pow(10, 15+i*0.1) * eV;
		candidate.current.setEnergy(E);
		epp2.process(&candidate, secondaries);
		EXPECT_TRUE(candidate.current.getEnergy() < E);
	}

	ElectronPairProduction epp3(ElectronPairProduction::CMBIR);
	for (int i = 0; i<80; i++) {
		double E = pow(10, 15+i*0.1) * eV;
		candidate.current.setEnergy(E);
		epp3.process(&candidate, secondaries);
		EXPECT_TRUE(candidate.current.getEnergy() < E);
	}
}

TEST(testElectronPairProduction, BelowEnergyTreshold) {
	// test if nothing happens below 1e15 eV
	ElectronPairProduction epp(ElectronPairProduction::CMB);

	Candidate candidate;
	candidate.current.setId(1000010010); // proton
	std::vector<Candidate *> secondaries;

	double E = 1e14 * eV;
	candidate.current.setEnergy(E);
	epp.process(&candidate, secondaries);
	EXPECT_DOUBLE_EQ(candidate.current.getEnergy(), E);
}

TEST(testElectronPairProduction, EnergyLossValues) {
	// test if energy loss corresponds to table
	std::vector<double> x;
	std::vector<double> y;
	std::ifstream infile("data/epair_cmbir.table");
	char header[256];
	for (int i = 0; i < 3; i++)
		infile.getline(header, 255);
	double a, b;
	while (infile.good()) {
		infile >> a >> b;
		x.push_back(a * eV);
		y.push_back(b * eV / Mpc);
	}
	infile.close();

	Candidate candidate;
	candidate.setCurrentStep(1 * Mpc);
	candidate.current.setId(1000010010); // proton
	std::vector<Candidate *> secondaries;

	ElectronPairProduction epp;
	for (int i = 0; i<x.size(); i++) {
		candidate.current.setEnergy(x[i]);
		epp.process(&candidate, secondaries);
		double dE = x[i] - candidate.current.getEnergy();
		double dE_table = y[i] * 1 * Mpc;
		EXPECT_NEAR(dE, dE_table, 1e-12);
	}
}

//TEST(testDecay, NeutronDecay) {
//	// test neutron life time
//	Candidate candidate;
//	candidate.setCurrentStep(1 * Mpc);
//	candidate.current.setId(1000010010); // proton
//	std::vector<Candidate *> secondaries;
//
//
//}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace mpc

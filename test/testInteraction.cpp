#include "mpc/Candidate.h"
#include "mpc/module/ElectronPairProduction.h"
#include "mpc/module/NuclearDecay.h"
#include "mpc/module/PhotoDisintegration.h"
#include "mpc/module/PhotoPionProduction.h"

#include "gtest/gtest.h"
#include <fstream>

namespace mpc {

TEST(ElectronPairProduction, EnergyDecreasing) {
	// test if energy loss occurs for protons with energies from 1e15 - 1e23 eV
	Candidate candidate;
	candidate.setCurrentStep(1 * Mpc);
	candidate.current.setId(getNucleusId(1, 1)); // proton

	ElectronPairProduction epp1(CMB);
	for (int i = 0; i < 80; i++) {
		double E = pow(10, 15 + i * 0.1) * eV;
		candidate.current.setEnergy(E);
		epp1.process(&candidate);
		EXPECT_TRUE(candidate.current.getEnergy() <= E);
	}

	ElectronPairProduction epp2(IRB);
	for (int i = 0; i < 80; i++) {
		double E = pow(10, 15 + i * 0.1) * eV;
		candidate.current.setEnergy(E);
		epp2.process(&candidate);
		EXPECT_TRUE(candidate.current.getEnergy() < E);
	}

	ElectronPairProduction epp3(CMB_IRB);
	for (int i = 0; i < 80; i++) {
		double E = pow(10, 15 + i * 0.1) * eV;
		candidate.current.setEnergy(E);
		epp3.process(&candidate);
		EXPECT_TRUE(candidate.current.getEnergy() < E);
	}
}

TEST(ElectronPairProduction, BelowEnergyTreshold) {
	// test if nothing happens below 1e15 eV
	ElectronPairProduction epp(CMB);
	Candidate candidate;
	candidate.current.setId(getNucleusId(1, 1)); // proton
	double E = 1e14 * eV;
	candidate.current.setEnergy(E);
	epp.process(&candidate);
	EXPECT_DOUBLE_EQ(candidate.current.getEnergy(), E);
}

TEST(ElectronPairProduction, EnergyLossValues) {
	// test if energy loss corresponds to table
	std::vector<double> x;
	std::vector<double> y;
	std::ifstream infile("data/ElectronPairProduction/cmbir.txt");
	while (infile.good()) {
		if (infile.peek() != '#') {
			double a, b;
			infile >> a >> b;
			if (infile) {
				x.push_back(a * eV);
				y.push_back(b * eV / Mpc);
			}
		}
		infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
	infile.close();

	Candidate candidate;
	candidate.setCurrentStep(1 * Mpc);
	candidate.current.setId(getNucleusId(1, 1)); // proton

	ElectronPairProduction epp;
	for (int i = 0; i < x.size(); i++) {
		candidate.current.setEnergy(x[i]);
		epp.process(&candidate);
		double dE = x[i] - candidate.current.getEnergy();
		double dE_table = y[i] * 1 * Mpc;
		EXPECT_NEAR(dE, dE_table, 1e-12);
	}
}



TEST(NuclearDecay, Neutron) {
	// test beta- decay n -> p
	Candidate candidate;
	candidate.setCurrentStep(1 * Mpc);
	candidate.current.setId(getNucleusId(1, 0));
	candidate.current.setEnergy(1 * EeV);
	NuclearDecay d;
	d.process(&candidate);
	EXPECT_EQ(candidate.current.getId(), getNucleusId(1,1));
}

TEST(NuclearDecay, Scandium44) {
	// test beta+ decay 44Sc -> 44Ca
	Candidate candidate;
	candidate.setCurrentStep(1 * Mpc);
	candidate.current.setId(getNucleusId(44, 21));
	candidate.current.setEnergy(1 * EeV);
	NuclearDecay d;
	d.process(&candidate);
	EXPECT_EQ(candidate.current.getId(), getNucleusId(44,20));
}

TEST(NuclearDecay, Li4) {
	// test proton dripping Li-4 -> He-3
	Candidate candidate;
	candidate.setCurrentStep(1 * kpc);
	candidate.current.setId(getNucleusId(4, 3));
	candidate.current.setEnergy(1 * EeV);
	NuclearDecay d;
	d.process(&candidate);
	EXPECT_EQ(getNucleusId(3,2), candidate.current.getId());
}

TEST(NuclearDecay, He5) {
	// test neturon dripping He-5 -> He-4
	Candidate candidate;
	candidate.setCurrentStep(1 * Mpc);
	candidate.current.setId(getNucleusId(5, 2));
	candidate.current.setEnergy(1 * EeV);
	NuclearDecay d;
	d.process(&candidate);
	EXPECT_EQ(getNucleusId(4,2), candidate.current.getId());
}

TEST(NuclearDecay, LimitNextStep) {
	// test if nextStep is limited
	// the decay length of a neutron at 10 EeV is 93 kpc
	NuclearDecay d;
	Candidate candidate;
	candidate.setNextStep(10 * Mpc);
	candidate.current.setId(getNucleusId(1, 0));
	candidate.current.setEnergy(10 * EeV);
	d.process(&candidate);
	EXPECT_TRUE(candidate.getNextStep() < 10 * Mpc);
}



TEST(PhotoDisintegration, Backgrounds) {
	PhotoPionProduction ppp1(CMB);
	PhotoPionProduction ppp2(IRB);
	PhotoPionProduction ppp3(CMB_IRB);
}

TEST(PhotoDisintegration, Carbon) {
	PhotoDisintegration pd;
	Candidate candidate;
	candidate.current.setId(getNucleusId(12, 6));
	candidate.current.setEnergy(200 * EeV);
	candidate.setCurrentStep(50 * Mpc);
	pd.process(&candidate);
	EXPECT_TRUE(candidate.current.getMassNumber() < 12);
	EXPECT_TRUE(candidate.current.getEnergy() < 200 * EeV);
	EXPECT_TRUE(candidate.secondaries.size() > 0);
}

TEST(PhotoDisintegration, Iron) {
	PhotoDisintegration pd;
	Candidate candidate;
	candidate.current.setId(getNucleusId(56, 26));
	candidate.current.setEnergy(200 * EeV);
	candidate.setCurrentStep(50 * Mpc);
	pd.process(&candidate);
	EXPECT_TRUE(candidate.current.getMassNumber() < 56);
	EXPECT_TRUE(candidate.current.getEnergy() < 200 * EeV);
	EXPECT_TRUE(candidate.secondaries.size() > 0);
}

TEST(PhotoDisintegration, LimitNextStep) {
	// test if nextStep is limited
	PhotoDisintegration pd;
	Candidate candidate;
	candidate.setNextStep(100 * Mpc);
	candidate.current.setId(getNucleusId(4, 2));
	candidate.current.setEnergy(200 * EeV);
	pd.process(&candidate);
	EXPECT_TRUE(candidate.getNextStep() < 100 * Mpc);
}



TEST(PhotoPionProduction, Backgrounds) {
	PhotoPionProduction ppp1(CMB);
	PhotoPionProduction ppp2(IRB);
	PhotoPionProduction ppp3(CMB_IRB);
}

TEST(PhotoPionProduction, Proton) {
	PhotoPionProduction ppp;
	Candidate candidate;
	candidate.setCurrentStep(100 * Mpc);
	candidate.current.setId(getNucleusId(1, 1));
	candidate.current.setEnergy(100 * EeV);
	ppp.process(&candidate);
	EXPECT_TRUE(candidate.current.getEnergy() / EeV < 100);
	EXPECT_EQ(candidate.current.getMassNumber(), 1);
}

TEST(PhotoPionProduction, Helium) {
	PhotoPionProduction ppp;
	Candidate candidate;
	candidate.setCurrentStep(100 * Mpc);
	candidate.current.setId(getNucleusId(4, 2));
	candidate.current.setEnergy(400 * EeV);
	ppp.process(&candidate);
	EXPECT_TRUE(candidate.current.getEnergy() / EeV < 300);
	EXPECT_TRUE(candidate.current.getMassNumber() < 4);
	EXPECT_TRUE(candidate.secondaries.size() > 0);
}

TEST(PhotoPionProduction, LimitNextStep) {
	// test if nextStep is limited
	PhotoPionProduction ppp;
	Candidate candidate;
	candidate.setNextStep(100 * Mpc);
	candidate.current.setId(getNucleusId(1, 1));
	candidate.current.setEnergy(200 * EeV);
	ppp.process(&candidate);
	EXPECT_TRUE(candidate.getNextStep() < 100 * Mpc);
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace mpc

#include "gtest/gtest.h"
#include "mpc/Candidate.h"
#include "mpc/module/ElectronPairProduction.h"
#include "mpc/module/NuclearDecay.h"
#include "mpc/module/PhotoDisintegration.h"
#include "mpc/module/PhotoPionProduction.h"

#include <fstream>

namespace mpc {

class ModuleTest: public testing::Test {
protected:
	Candidate candidate;
	std::vector<Candidate *> secondaries;
};

TEST_F(ModuleTest, ElectronPairProduction_EnergyDecreasing) {
	// test if energy loss occurs for protons with energies from 1e15 - 1e23 eV
	candidate.setCurrentStep(1 * Mpc);
	candidate.current.setId(getNucleusId(1, 1)); // proton

	ElectronPairProduction epp1(ElectronPairProduction::CMB);
	for (int i = 0; i < 80; i++) {
		double E = pow(10, 15 + i * 0.1) * eV;
		candidate.current.setEnergy(E);
		epp1.process(&candidate, secondaries);
		EXPECT_TRUE(candidate.current.getEnergy() <= E);
	}

	ElectronPairProduction epp2(ElectronPairProduction::IR);
	for (int i = 0; i < 80; i++) {
		double E = pow(10, 15 + i * 0.1) * eV;
		candidate.current.setEnergy(E);
		epp2.process(&candidate, secondaries);
		EXPECT_TRUE(candidate.current.getEnergy() < E);
	}

	ElectronPairProduction epp3(ElectronPairProduction::CMBIR);
	for (int i = 0; i < 80; i++) {
		double E = pow(10, 15 + i * 0.1) * eV;
		candidate.current.setEnergy(E);
		epp3.process(&candidate, secondaries);
		EXPECT_TRUE(candidate.current.getEnergy() < E);
	}
}

TEST_F(ModuleTest, ElectronPairProduction_BelowEnergyTreshold) {
	// test if nothing happens below 1e15 eV
	ElectronPairProduction epp(ElectronPairProduction::CMB);
	candidate.current.setId(getNucleusId(1, 1)); // proton
	double E = 1e14 * eV;
	candidate.current.setEnergy(E);
	epp.process(&candidate, secondaries);
	EXPECT_DOUBLE_EQ(candidate.current.getEnergy(), E);
}

TEST_F(ModuleTest, ElectronPairProduction_EnergyLossValues) {
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

	candidate.setCurrentStep(1 * Mpc);
	candidate.current.setId(getNucleusId(1, 1)); // proton

	ElectronPairProduction epp;
	for (int i = 0; i < x.size(); i++) {
		candidate.current.setEnergy(x[i]);
		epp.process(&candidate, secondaries);
		double dE = x[i] - candidate.current.getEnergy();
		double dE_table = y[i] * 1 * Mpc;
		EXPECT_NEAR(dE, dE_table, 1e-12);
	}
}

TEST_F(ModuleTest, NuclearDecay_Neutron) {
	// test beta- decay n -> p
	candidate.setCurrentStep(1 * Mpc);
	candidate.current.setId(getNucleusId(1, 0));
	candidate.current.setEnergy(1 * EeV);
	NuclearDecay d;
	d.process(&candidate, secondaries);
	EXPECT_EQ(candidate.current.getId(), getNucleusId(1,1));
}

TEST_F(ModuleTest, NuclearDecay_Scandium44) {
	// test beta+ decay 44Sc -> 44Ca
	candidate.setCurrentStep(1 * Mpc);
	candidate.current.setId(getNucleusId(44, 21));
	candidate.current.setEnergy(1 * EeV);
	NuclearDecay d;
	d.process(&candidate, secondaries);
	EXPECT_EQ(candidate.current.getId(), getNucleusId(44,20));
}

TEST_F(ModuleTest, NuclearDecay_LimitNextStep) {
	NuclearDecay d;
	candidate.setNextStep(10 * Mpc);
	candidate.current.setId(getNucleusId(1, 0));
	candidate.current.setEnergy(10 * EeV);
	d.process(&candidate, secondaries);
	EXPECT_TRUE(candidate.getNextStep() < 10 * Mpc);
}



TEST_F(ModuleTest, PhotoDisintegration_Carbon) {
	PhotoDisintegration pd;
	candidate.current.setId(getNucleusId(12, 6));
	candidate.current.setEnergy(200 * EeV);
	candidate.setCurrentStep(50 * Mpc);
	pd.process(&candidate, secondaries);
	EXPECT_TRUE(candidate.current.getMassNumber() < 12);
	EXPECT_TRUE(candidate.current.getEnergy() < 200 * EeV);
}

TEST_F(ModuleTest, PhotoDisintegration_Iron) {
	PhotoDisintegration pd;
	candidate.current.setId(getNucleusId(56, 26));
	candidate.current.setEnergy(200 * EeV);
	candidate.setCurrentStep(50 * Mpc);
	pd.process(&candidate, secondaries);
	EXPECT_TRUE(candidate.current.getMassNumber() < 56);
	EXPECT_TRUE(candidate.current.getEnergy() < 200 * EeV);
}



TEST_F(ModuleTest, PhotoPionProduction_Proton) {
	PhotoPionProduction ppp(PhotoPionProduction::CMBIR);
	candidate.setCurrentStep(100 * Mpc);
	candidate.current.setId(getNucleusId(1, 1));
	candidate.current.setEnergy(100 * EeV);
	ppp.process(&candidate, secondaries);
	EXPECT_TRUE(candidate.current.getEnergy() / EeV < 100);
	EXPECT_EQ(candidate.current.getMassNumber(), 1);
}

TEST_F(ModuleTest, PhotoPionProduction_Iron) {
	PhotoPionProduction ppp(PhotoPionProduction::CMBIR);
	candidate.setCurrentStep(100 * Mpc);
	candidate.current.setId(getNucleusId(56, 26));
	candidate.current.setEnergy(100 * EeV);
	ppp.process(&candidate, secondaries);
	EXPECT_TRUE(candidate.current.getEnergy() / EeV < 1000);
	EXPECT_TRUE(candidate.current.getMassNumber() == 56);
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace mpc

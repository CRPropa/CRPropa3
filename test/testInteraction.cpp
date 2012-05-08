#include "mpc/Candidate.h"
#include "mpc/module/ElectronPairProduction.h"
#include "mpc/module/NuclearDecay.h"
#include "mpc/module/PhotoDisintegration.h"
#include "mpc/module/PhotoPionProduction.h"

#include "gtest/gtest.h"
#include <fstream>

namespace mpc {

TEST(ElectronPairProduction, EnergyDecreasing) {
	// Test if energy loss occurs for protons with energies from 1e15 - 1e23 eV
	Candidate c;
	c.setCurrentStep(1 * Mpc);
	c.current.setId(getNucleusId(1, 1)); // proton

	ElectronPairProduction epp1(CMB);
	for (int i = 0; i < 80; i++) {
		double E = pow(10, 15 + i * 0.1) * eV;
		c.current.setEnergy(E);
		epp1.process(&c);
		EXPECT_TRUE(c.current.getEnergy() <= E);
	}

	ElectronPairProduction epp2(IRB);
	for (int i = 0; i < 80; i++) {
		double E = pow(10, 15 + i * 0.1) * eV;
		c.current.setEnergy(E);
		epp2.process(&c);
		EXPECT_TRUE(c.current.getEnergy() < E);
	}

	ElectronPairProduction epp3(CMB_IRB);
	for (int i = 0; i < 80; i++) {
		double E = pow(10, 15 + i * 0.1) * eV;
		c.current.setEnergy(E);
		epp3.process(&c);
		EXPECT_TRUE(c.current.getEnergy() < E);
	}
}

TEST(ElectronPairProduction, BelowEnergyTreshold) {
	// Test if nothing happens below 1e15 eV
	ElectronPairProduction epp(CMB);
	Candidate c;
	c.current.setId(getNucleusId(1, 1)); // proton
	double E = 1e14 * eV;
	c.current.setEnergy(E);
	epp.process(&c);
	EXPECT_DOUBLE_EQ(c.current.getEnergy(), E);
}

TEST(ElectronPairProduction, EnergyLossValues) {
	// Test if energy loss corresponds to the data table.
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

	Candidate c;
	c.setCurrentStep(1 * Mpc);
	c.current.setId(getNucleusId(1, 1)); // proton

	ElectronPairProduction epp;
	for (int i = 0; i < x.size(); i++) {
		c.current.setEnergy(x[i]);
		epp.process(&c);
		double dE = x[i] - c.current.getEnergy();
		double dE_table = y[i] * 1 * Mpc;
		EXPECT_NEAR(dE, dE_table, 1e-12);
	}
}

TEST(NuclearDecay, Neutron) {
	// Quantitative test of decaying neutrons at rest.
	// The mean decay time is expected to be within 1% of the literature value.
	// This test can stochastically fail.
	Candidate candidate;
	candidate.current.setId(getNucleusId(1, 0));
	candidate.current.setEnergy(mass_neutron * c_squared);
	NuclearDecay decay;
	InteractionState state;
	double tau = 0;
	for (int i = 0; i < 100000; i++) {
		decay.setNextInteraction(&candidate, state);
		tau += state.distance;
	}
	tau /= c_light * 100000;
	EXPECT_NEAR(tau, 881.46, 8.8);
}

TEST(NuclearDecay, Scandium44) {
	// Test beta+ decay of 44Sc 44Ca.
	// This test can stochastically fail.
	NuclearDecay d(true, true);
	Candidate c;
	c.setCurrentStep(1 * Mpc);
	c.current.setId(getNucleusId(44, 21));
	c.current.setEnergy(1 * EeV);
	double gamma = c.current.getLorentzFactor();
	d.process(&c);
	// primary
	EXPECT_EQ(getNucleusId(44,20), c.current.getId());
	EXPECT_DOUBLE_EQ(gamma, c.current.getLorentzFactor());
	// secondaries
	EXPECT_EQ(2, c.secondaries.size());
	Candidate c1 = *c.secondaries[0];
	Candidate c2 = *c.secondaries[1];
	EXPECT_EQ(-11, c1.current.getId());
	// positron
	EXPECT_EQ( 12, c2.current.getId());
	// electron neutrino
}

TEST(NuclearDecay, Li4) {
	// Test proton dripping of Li-4 to He-3.
	// This test can stochastically fail.
	NuclearDecay d;
	Candidate c;
	c.setCurrentStep(1 * kpc);
	c.current.setId(getNucleusId(4, 3));
	c.current.setEnergy(4 * EeV);
	d.process(&c);
	// primary
	EXPECT_EQ(getNucleusId(3,2), c.current.getId());
	EXPECT_EQ(1, c.secondaries.size());
	// secondary
	Candidate c1 = *c.secondaries[0];
	EXPECT_EQ(getNucleusId(1,1), c1.current.getId());
	EXPECT_EQ(1, c1.current.getEnergy() / EeV);
}

TEST(NuclearDecay, He5) {
	// Test neturon dripping of He-5 to He-4.
	// This test can stochastically fail if no interaction occurs over 1 Mpc.
	NuclearDecay d;
	Candidate c;
	c.setCurrentStep(1 * Mpc);
	c.current.setId(getNucleusId(5, 2));
	c.current.setEnergy(5 * EeV);
	d.process(&c);
	// primary
	EXPECT_EQ(getNucleusId(4,2), c.current.getId());
	EXPECT_EQ(4, c.current.getEnergy() / EeV);
	// secondary
	Candidate c2 = *c.secondaries[0];
	EXPECT_EQ(getNucleusId(1,0), c2.current.getId());
	EXPECT_EQ(1, c2.current.getEnergy() / EeV);
}

TEST(NuclearDecay, LimitNextStep) {
	// Test if next step is limited.
	NuclearDecay d;
	Candidate c;
	c.setNextStep(std::numeric_limits<double>::max());
	c.current.setId(getNucleusId(1, 0));
	c.current.setEnergy(10 * EeV);
	d.process(&c);
	EXPECT_LT(c.getNextStep(), std::numeric_limits<double>::max());
}

TEST(NuclearDecay, AllWorking) {
	// Test if all nuclear decays are working.
	NuclearDecay d;
	Candidate c;
	InteractionState interaction;

	std::ifstream infile(getDataPath("/NuclearDecay/decayTable.txt").c_str());
	while (infile.good()) {
		if (infile.peek() != '#') {
			int Z, N, channel, foo;
			infile >> Z >> N >> interaction.channel >> foo;

			c.current.setId(getNucleusId(Z + N, Z));
			c.current.setEnergy(80 * EeV);
			c.setInteractionState(d.getDescription(), interaction);
			d.performInteraction(&c);
		}
		infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
	infile.close();
}

TEST(PhotoDisintegration, Carbon) {
	// Test if a 100 EeV C-12 nucleus photo-disintegrates (at least once) over a distance of 50 Mpc.
	// This test can stochastically fail if no interaction occurs over 50 Mpc.
	PhotoDisintegration pd;
	Candidate c;
	c.current.setId(getNucleusId(12, 6));
	c.current.setEnergy(100 * EeV);
	c.setCurrentStep(50 * Mpc);
	pd.process(&c);

	EXPECT_TRUE(c.current.getEnergy() < 100 * EeV);
	// energy loss
	EXPECT_TRUE(c.secondaries.size() > 0);
	// secondaries produced

	int A = c.current.getMassNumber();
	int Z = c.current.getChargeNumber();
	double E = c.current.getEnergy();

	for (int i = 0; i < c.secondaries.size(); i++) {
		A += (*c.secondaries[i]).current.getMassNumber();
		Z += (*c.secondaries[i]).current.getChargeNumber();
		E += (*c.secondaries[i]).current.getEnergy();
	}
	EXPECT_EQ(12, A);
	// nucleon number conserved
	EXPECT_EQ(6, Z);
	// proton number conserved
	EXPECT_DOUBLE_EQ(100 * EeV, E);
	// energy conserved
}

TEST(PhotoDisintegration, Iron) {
	// Test if a 100 EeV Fe-56 nucleus photo-disintegrates (at least once) over a distance of 50 Mpc.
	// This test can stochastically fail if no interaction occurs over 50 Mpc.
	PhotoDisintegration pd(IRB);
	Candidate c;
	c.current.setId(getNucleusId(56, 26));
	c.current.setEnergy(100 * EeV);
	c.setCurrentStep(50 * Mpc);
	pd.process(&c);

	EXPECT_TRUE(c.current.getEnergy() < 100 * EeV);
	// energy loss
	EXPECT_TRUE(c.secondaries.size() > 0);
	// secondaries produced

	int A = c.current.getMassNumber();
	int Z = c.current.getChargeNumber();
	double E = c.current.getEnergy();

	for (int i = 0; i < c.secondaries.size(); i++) {
		A += (*c.secondaries[i]).current.getMassNumber();
		Z += (*c.secondaries[i]).current.getChargeNumber();
		E += (*c.secondaries[i]).current.getEnergy();
	}
	EXPECT_EQ(56, A);
	// nucleon number conserved
	EXPECT_EQ(26, Z);
	// proton number conserved
	EXPECT_DOUBLE_EQ(100 * EeV, E);
	// energy conserved
}

TEST(PhotoDisintegration, LimitNextStep) {
	// Test if the interaction limits the next propagation step.
	PhotoDisintegration pd;
	Candidate c;
	c.setNextStep(std::numeric_limits<double>::max());
	c.current.setId(getNucleusId(4, 2));
	c.current.setEnergy(200 * EeV);
	pd.process(&c);
	EXPECT_LT(c.getNextStep(), std::numeric_limits<double>::max());
}

TEST(PhotoDisintegration, AllWorkingCMB) {
	// Test if all photo-disintegrations are working.
	PhotoDisintegration pd(CMB);
	Candidate c;
	InteractionState interaction;

	std::ifstream infile(
			getDataPath("PhotoDisintegration/PDtable_CMB.txt").c_str());
	std::string line;
	while (std::getline(infile, line)) {
		if (line[0] == '#')
			continue;
		std::stringstream lineStream(line);
		int Z, N;
		lineStream >> Z;
		lineStream >> N;
		lineStream >> interaction.channel;

		double y;
		for (size_t i = 0; i < 200; i++) {
			lineStream >> y;
			EXPECT_TRUE(lineStream);
			// test if all 200 entries are present
		}

		c.current.setId(getNucleusId(Z + N, Z));
		c.current.setEnergy(80 * EeV);
		c.setInteractionState(pd.getDescription(), interaction);
		pd.performInteraction(&c);
	}
	infile.close();
}

TEST(PhotoDisintegration, AllWorkingIRB) {
	// Test if all photo-disintegrations are working.
	PhotoDisintegration pd(IRB);
	Candidate c;
	InteractionState interaction;

	std::ifstream infile(
			getDataPath("PhotoDisintegration/PDtable_IRB.txt").c_str());
	std::string line;
	while (std::getline(infile, line)) {
		if (line[0] == '#')
			continue;
		std::stringstream lineStream(line);
		int Z, N;
		lineStream >> Z;
		lineStream >> N;
		lineStream >> interaction.channel;

		double y;
		for (size_t i = 0; i < 200; i++) {
			lineStream >> y;
			EXPECT_TRUE(lineStream);
			// test if all 200 entries are present
		}

		c.current.setId(getNucleusId(Z + N, Z));
		c.current.setEnergy(80 * EeV);
		c.setInteractionState(pd.getDescription(), interaction);
		pd.performInteraction(&c);
	}
	infile.close();
}

TEST(PhotoPionProduction, Backgrounds) {
	// Test if interaction data files are loaded.
	PhotoPionProduction ppp1(CMB);
	PhotoPionProduction ppp2(IRB);
}

TEST(PhotoPionProduction, Proton) {
	// Test photo-pion interaction for 100 EeV proton.
	// This test can stochastically fail if no interaction occurs over 100 Mpc.
	PhotoPionProduction ppp;
	Candidate c;
	c.setCurrentStep(100 * Mpc);
	c.current.setId(getNucleusId(1, 1));
	c.current.setEnergy(100 * EeV);
	ppp.process(&c);
	EXPECT_TRUE(c.current.getEnergy() / EeV < 100);
	// energy loss
	EXPECT_EQ(1, c.current.getMassNumber());
	// nucleon number conserved
	EXPECT_EQ(0, c.secondaries.size());
	// no (nucleonic) secondaries
}

TEST(PhotoPionProduction, Helium) {
	// Test photo-pion interaction for 400 EeV He nucleus.
	// This test can stochastically fail if no interaction occurs over 100 Mpc.
	PhotoPionProduction ppp;
	Candidate c;
	c.setCurrentStep(100 * Mpc);
	c.current.setId(getNucleusId(4, 2));
	c.current.setEnergy(400 * EeV);
	ppp.process(&c);
	EXPECT_LT(c.current.getEnergy(), 400 * EeV);
	EXPECT_TRUE(c.current.getMassNumber() < 4);
	EXPECT_TRUE(c.secondaries.size() > 0);
}

TEST(PhotoPionProduction, LimitNextStep) {
	// Test if the interaction limits the next propagation step.
	PhotoPionProduction ppp;
	Candidate c;
	c.setNextStep(std::numeric_limits<double>::max());
	c.current.setId(getNucleusId(1, 1));
	c.current.setEnergy(200 * EeV);
	ppp.process(&c);
	EXPECT_LT(c.getNextStep(), std::numeric_limits<double>::max());
}

TEST(SophiaPhotoPionProduction, withoutSecondaries) {
	// Test photo-pion (SOPHIA) interaction for 100 EeV proton.
	// This test can stochastically fail if no interaction occurs over 100 Mpc.
	SophiaPhotoPionProduction ppp;
	Candidate c;
	c.setCurrentStep(100 * Mpc);
	c.current.setId(getNucleusId(1, 1));
	c.current.setEnergy(100 * EeV);
	ppp.process(&c);
	EXPECT_GT(100 * EeV, c.current.getEnergy());
	// energy loss
	EXPECT_EQ(1, c.current.getMassNumber());
	// nucleon number conserved
	EXPECT_EQ(0, c.secondaries.size());
	// secondaries turned off
}

TEST(SophiaPhotoPionProduction, withSecondaries) {
	// Test photo-pion interaction for 100 EeV proton.
	// This test can stochastically fail if no interaction occurs over 100 Mpc.
	SophiaPhotoPionProduction ppp(CMB, true, true, true);
	Candidate c;
	c.current.setId(getNucleusId(1, 1));
	c.current.setEnergy(100 * EeV);
	InteractionState interaction;
	ppp.setNextInteraction(&c, interaction);
	ppp.performInteraction(&c);
	EXPECT_GT(c.secondaries.size(), 1);
	// secondaries turned on
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace mpc

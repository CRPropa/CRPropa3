#include "crpropa/Candidate.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include <HepPID/ParticleIDMethods.hh>
#include "crpropa/module/ElectronPairProduction.h"
#include "crpropa/module/NuclearDecay.h"
#include "crpropa/module/PhotoDisintegration.h"
#include "crpropa/module/PhotoPionProduction.h"
#include "crpropa/module/Redshift.h"
#include "gtest/gtest.h"

#include <fstream>

namespace crpropa {

TEST(ElectronPairProduction, energyDecreasing) {
	// Test if energy loss occurs for protons with energies from 1e15 - 1e23 eV
	Candidate c;
	c.setCurrentStep(2 * Mpc);
	c.current.setId(nucleusId(1, 1)); // proton

	ElectronPairProduction epp1(CMB);
	for (int i = 0; i < 80; i++) {
		double E = pow(10, 15 + i * 0.1) * eV;
		c.current.setEnergy(E);
		epp1.process(&c);
		EXPECT_LE(c.current.getEnergy(), E);
	}

	ElectronPairProduction epp2(IRB);
	for (int i = 0; i < 80; i++) {
		double E = pow(10, 15 + i * 0.1) * eV;
		c.current.setEnergy(E);
		epp2.process(&c);
		EXPECT_LE(c.current.getEnergy(), E);
	}
}

TEST(ElectronPairProduction, belowEnergyTreshold) {
	// Test if nothing happens below 1e15 eV
	ElectronPairProduction epp(CMB);
	Candidate c;
	c.current.setId(nucleusId(1, 1)); // proton
	double E = 1e14 * eV;
	c.current.setEnergy(E);
	epp.process(&c);
	EXPECT_DOUBLE_EQ(E, c.current.getEnergy());
}

TEST(ElectronPairProduction, thisIsNotNucleonic) {
	// Test if non-nuclei are skipped
	ElectronPairProduction epp(CMB);
	Candidate c;
	c.current.setId(11); // electron
	double E = 1e20 * eV;
	c.current.setEnergy(E);
	epp.process(&c);
	EXPECT_DOUBLE_EQ(E, c.current.getEnergy());
}

TEST(ElectronPairProduction, valuesCMB) {
	// Test if energy loss corresponds to the data table.
	std::vector<double> x;
	std::vector<double> y;
	std::ifstream infile(getDataPath("pair_CMB.txt").c_str());
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
	c.current.setId(nucleusId(1, 1)); // proton

	ElectronPairProduction epp(CMB);
	for (int i = 0; i < x.size(); i++) {
		c.current.setEnergy(x[i]);
		epp.process(&c);
		double dE = x[i] - c.current.getEnergy();
		double dE_table = y[i] * 1 * Mpc;
		EXPECT_NEAR(dE_table, dE, 1e-12);
	}
}

TEST(ElectronPairProduction, valuesIRB) {
	// Test if energy loss corresponds to the data table.
	std::vector<double> x;
	std::vector<double> y;
	std::ifstream infile(getDataPath("pairIRB.txt").c_str());
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
	c.current.setId(nucleusId(1, 1)); // proton

	ElectronPairProduction epp(IRB);
	for (int i = 0; i < x.size(); i++) {
		c.current.setEnergy(x[i]);
		epp.process(&c);
		double dE = x[i] - c.current.getEnergy();
		double dE_table = y[i] * 1 * Mpc;
		EXPECT_NEAR(dE, dE_table, 1e-12);
	}
}

TEST(NuclearDecay, scandium44) {
	// Test beta+ decay of 44Sc 44Ca.
	// This test can stochastically fail.
	NuclearDecay d(true, true);
	Candidate c;
	c.setCurrentStep(1 * Mpc);
	c.current.setId(nucleusId(44, 21));
	c.current.setEnergy(1 * EeV);
	double gamma = c.current.getLorentzFactor();
	d.process(&c);
	// primary
	EXPECT_EQ(nucleusId(44, 20), c.current.getId());
	EXPECT_DOUBLE_EQ(gamma, c.current.getLorentzFactor());
	// secondaries
	EXPECT_EQ(2, c.secondaries.size());
	Candidate c1 = *c.secondaries[0];
	Candidate c2 = *c.secondaries[1];
	EXPECT_EQ(-11, c1.current.getId());
	// positron
	EXPECT_EQ(12, c2.current.getId());
	// electron neutrino
}

TEST(NuclearDecay, lithium4) {
	// Test proton dripping of Li-4 to He-3.
	// This test can stochastically fail.
	NuclearDecay d;
	Candidate c;
	c.setCurrentStep(1 * kpc);
	c.current.setId(nucleusId(4, 3));
	c.current.setEnergy(4 * EeV);
	d.process(&c);
	// primary
	EXPECT_EQ(nucleusId(3, 2), c.current.getId());
	EXPECT_EQ(1, c.secondaries.size());
	// secondary
	Candidate c1 = *c.secondaries[0];
	EXPECT_EQ(nucleusId(1, 1), c1.current.getId());
	EXPECT_EQ(1, c1.current.getEnergy() / EeV);
}

TEST(NuclearDecay, helium5) {
	// Test neturon dripping of He-5 to He-4.
	// This test can stochastically fail if no interaction occurs over 1 Mpc.
	NuclearDecay d;
	Candidate c;
	c.setCurrentStep(1 * Mpc);
	c.current.setId(nucleusId(5, 2));
	c.current.setEnergy(5 * EeV);
	d.process(&c);
	// primary
	EXPECT_EQ(nucleusId(4, 2), c.current.getId());
	EXPECT_EQ(4, c.current.getEnergy() / EeV);
	// secondary
	Candidate c2 = *c.secondaries[0];
	EXPECT_EQ(nucleusId(1, 0), c2.current.getId());
	EXPECT_EQ(1, c2.current.getEnergy() / EeV);
}

TEST(NuclearDecay, limitNextStep) {
	// Test if next step is limited in case of a neutron.
	NuclearDecay decay;
	Candidate c;
	c.setNextStep(std::numeric_limits<double>::max());
	c.current.setId(nucleusId(1, 0));
	c.current.setEnergy(10 * EeV);
	decay.process(&c);
	EXPECT_LT(c.getNextStep(), std::numeric_limits<double>::max());
}

TEST(NuclearDecay, allWorking) {
	// Test if all nuclear decays are working.
	NuclearDecay d;
	Candidate c;

	std::ifstream infile(getDataPath("nuclear_decay.txt").c_str());
	while (infile.good()) {
		if (infile.peek() != '#') {
			int Z, N, channel, foo;
			infile >> Z >> N >> channel >> foo;
			c.current.setId(nucleusId(Z + N, Z));
			c.current.setEnergy(80 * EeV);
			d.performInteraction(&c, channel);
		}
		infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
	infile.close();
}

TEST(NuclearDecay, thisIsNotNucleonic) {
	// Test if noting happens to an electron
	NuclearDecay decay;
	Candidate c;
	c.setNextStep(std::numeric_limits<double>::max());
	c.current.setId(11); // electron
	c.current.setEnergy(10 * EeV);
	decay.process(&c);
	EXPECT_EQ(11, c.current.getId());
	EXPECT_EQ(10 * EeV, c.current.getEnergy());
}

TEST(PhotoDisintegration, carbon) {
	// Test if a 100 EeV C-12 nucleus photo-disintegrates (at least once) over a distance of 50 Mpc.
	// This test can stochastically fail if no interaction occurs over 50 Mpc.
	PhotoDisintegration pd;
	Candidate c;
	int id = nucleusId(12, 6);
	c.current.setId(id);
	c.current.setEnergy(100 * EeV);
	c.setCurrentStep(50 * Mpc);
	pd.process(&c);

	EXPECT_TRUE(c.current.getEnergy() < 100 * EeV);
	// energy loss
	EXPECT_TRUE(c.secondaries.size() > 0);
	// secondaries produced

	double E = c.current.getEnergy();
	id = c.current.getId();
	int A = massNumber(id);
	int Z = chargeNumber(id);

	for (int i = 0; i < c.secondaries.size(); i++) {
		E += (*c.secondaries[i]).current.getEnergy();
		id = (*c.secondaries[i]).current.getId();
		A += massNumber(id);
		Z += chargeNumber(id);
	}
	EXPECT_EQ(12, A);
	// nucleon number conserved
	EXPECT_EQ(6, Z);
	// proton number conserved
	EXPECT_DOUBLE_EQ(100 * EeV, E);
	// energy conserved
}

TEST(PhotoDisintegration, iron) {
	// Test if a 100 EeV Fe-56 nucleus photo-disintegrates (at least once) over a distance of 50 Mpc.
	// This test can stochastically fail if no interaction occurs over 50 Mpc.
	PhotoDisintegration pd(IRB);
	Candidate c;
	int id = nucleusId(56, 26);
	c.current.setId(id);
	c.current.setEnergy(100 * EeV);
	c.setCurrentStep(50 * Mpc);
	pd.process(&c);

	EXPECT_TRUE(c.current.getEnergy() < 100 * EeV);
	// energy loss
	EXPECT_TRUE(c.secondaries.size() > 0);
	// secondaries produced

	double E = c.current.getEnergy();
	id = c.current.getId();
	int A = massNumber(id);
	int Z = chargeNumber(id);

	for (int i = 0; i < c.secondaries.size(); i++) {
		E += (*c.secondaries[i]).current.getEnergy();
		id = (*c.secondaries[i]).current.getId();
		A += massNumber(id);
		Z += chargeNumber(id);
	}
	EXPECT_EQ(56, A);
	// nucleon number conserved
	EXPECT_EQ(26, Z);
	// proton number conserved
	EXPECT_DOUBLE_EQ(100 * EeV, E);
	// energy conserved
}

TEST(PhotoDisintegration, thisIsNotNucleonic) {
	// Test that nothing happens to an electron
	PhotoDisintegration pd;
	Candidate c;
	c.setCurrentStep(1 * Mpc);
	c.current.setId(11); // electron
	c.current.setEnergy(10 * EeV);
	pd.process(&c);
	EXPECT_EQ(11, c.current.getId());
	EXPECT_EQ(10 * EeV, c.current.getEnergy());
}

TEST(PhotoDisintegration, limitNextStep) {
	// Test if the interaction limits the next propagation step.
	PhotoDisintegration pd;
	Candidate c;
	c.setNextStep(std::numeric_limits<double>::max());
	c.current.setId(nucleusId(4, 2));
	c.current.setEnergy(200 * EeV);
	pd.process(&c);
	EXPECT_LT(c.getNextStep(), std::numeric_limits<double>::max());
}

TEST(PhotoDisintegration, allIsotopes) {
	// Test if all isotopes are handled.
	PhotoDisintegration pd1(CMB);
	PhotoDisintegration pd2(IRB);
	Candidate c;
	c.setCurrentStep(10 * Mpc);

	for (int Z = 1; Z <= 26; Z++) {
		for (int N = 1; N <= 30; N++) {

			c.current.setId(nucleusId(Z + N, Z));
			c.current.setEnergy(80 * EeV);
			pd1.process(&c);

			c.current.setId(nucleusId(Z + N, Z));
			c.current.setEnergy(80 * EeV);
			pd2.process(&c);
		}
	}
}

TEST(PhotoPionProduction, backgrounds) {
	// Test if interaction data files are loaded.
	PhotoPionProduction ppp1(CMB);
	PhotoPionProduction ppp2(IRB);
}

TEST(PhotoPionProduction, proton) {
	// Test photo-pion interaction for 100 EeV proton.
	// This test can stochastically fail if no interaction occurs over 100 Mpc.
	PhotoPionProduction ppp;
	Candidate c;
	c.setCurrentStep(100 * Mpc);
	c.current.setId(nucleusId(1, 1));
	c.current.setEnergy(100 * EeV);
	ppp.process(&c);
	EXPECT_TRUE(c.current.getEnergy() / EeV < 100); // energy loss
	int id = c.current.getId();
	EXPECT_EQ(1, massNumber(id)); // nucleon number conserved
	EXPECT_EQ(0, c.secondaries.size()); // no (nucleonic) secondaries
}

TEST(PhotoPionProduction, helium) {
	// Test photo-pion interaction for 400 EeV He nucleus.
	// This test can stochastically fail if no interaction occurs over 100 Mpc.
	PhotoPionProduction ppp;
	Candidate c;
	c.setCurrentStep(100 * Mpc);
	c.current.setId(nucleusId(4, 2));
	c.current.setEnergy(400 * EeV);
	ppp.process(&c);
	EXPECT_LT(c.current.getEnergy(), 400 * EeV);
	int id = c.current.getId();
	EXPECT_TRUE(massNumber(id) < 4);
	EXPECT_TRUE(c.secondaries.size() > 0);
}

TEST(PhotoPionProduction, thisIsNotNucleonic) {
	// Test if noting happens to an electron
	PhotoPionProduction ppp;
	Candidate c;
	c.setCurrentStep(1 * Mpc);
	c.setNextStep(std::numeric_limits<double>::max());
	c.current.setId(11); // electron
	c.current.setEnergy(10 * EeV);
	ppp.process(&c);
	EXPECT_EQ(11, c.current.getId());
	EXPECT_EQ(10 * EeV, c.current.getEnergy());
}

TEST(PhotoPionProduction, limitNextStep) {
	// Test if the interaction limits the next propagation step.
	PhotoPionProduction ppp;
	Candidate c;
	c.setNextStep(std::numeric_limits<double>::max());
	c.current.setId(nucleusId(1, 1));
	c.current.setEnergy(200 * EeV);
	ppp.process(&c);
	EXPECT_LT(c.getNextStep(), std::numeric_limits<double>::max());
}

TEST(PhotoPionProduction, withoutSecondaries) {
	// Test photo-pion (SOPHIA) interaction for 100 EeV proton.
	// This test can stochastically fail if no interaction occurs over 100 Mpc.
	PhotoPionProduction ppp;
	Candidate c;
	c.setCurrentStep(100 * Mpc);
	c.current.setId(nucleusId(1, 1));
	c.current.setEnergy(100 * EeV);
	ppp.process(&c);

	// energy loss
	EXPECT_GT(100 * EeV, c.current.getEnergy());
	// nucleon number conserved
	int id = c.current.getId();
	EXPECT_EQ(1, massNumber(id));
	// secondaries turned off
	EXPECT_EQ(0, c.secondaries.size());
}

TEST(PhotoPionProduction, withSecondaries) {
	// Test photo-pion interaction for 100 EeV proton.
	// This test can stochastically fail if no interaction occurs over 100 Mpc.
	PhotoPionProduction ppp(CMB, true, true, true);
	Candidate c;
	c.current.setId(nucleusId(1, 1));
	c.current.setEnergy(100 * EeV);
	c.setCurrentStep(100 * Mpc);
	ppp.process(&c);
	// there should be secondaries
	EXPECT_GT(c.secondaries.size(), 1);
}

TEST(Redshift, simpleTest) {
	// Test if redshift is decreased and adiabatic energy loss is applied.
	Redshift redshift;

	Candidate c;
	c.setRedshift(0.024);
	c.current.setEnergy(100 * EeV);
	c.setCurrentStep(1 * Mpc);

	redshift.process(&c);
	EXPECT_GT(0.024, c.getRedshift()); // expect redshift decrease
	EXPECT_GT(100, c.current.getEnergy() / EeV); // expect energy loss
}

TEST(Redshift, limitRedshiftDecrease) {
	// Test if the redshift decrease is limited to z_updated >= 0.
	Redshift redshift;

	Candidate c;
	c.setRedshift(0.024); // roughly corresponds to 100 Mpc
	c.setCurrentStep(150 * Mpc);

	redshift.process(&c);
	EXPECT_DOUBLE_EQ(0, c.getRedshift());
}


TEST(PIDdigit, consistencyWithReferenceImplementation){
	// Tests the performance improved version against the default one
	unsigned long testPID = rand() % 1000000000 + 1000000000;
	for(size_t i=1; i < 8; i++)
	{
		HepPID::location loc = (HepPID::location) i;
		unsigned short newResult = HepPID::digit(loc, testPID);
		//original implementation
		int numerator = (int) std::pow(10.0,(loc-1));
		EXPECT_EQ(newResult, (HepPID::abspid(testPID)/numerator)%10);
	}
}


int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace crpropa

#include "crpropa/Candidate.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/module/ElectronPairProduction.h"
#include "crpropa/module/NuclearDecay.h"
#include "crpropa/module/PhotoDisintegration.h"
#include "crpropa/module/ElasticScattering.h"
#include "crpropa/module/PhotoPionProduction.h"
#include "crpropa/module/Redshift.h"
#include "crpropa/module/EMPairProduction.h"
#include "crpropa/module/EMDoublePairProduction.h"
#include "crpropa/module/EMTripletPairProduction.h"
#include "crpropa/module/EMInverseComptonScattering.h"
#include "gtest/gtest.h"

#include <fstream>

namespace crpropa {

// ElectronPairProduction -----------------------------------------------------
TEST(ElectronPairProduction, allBackgrounds) {
	// Test if interaction data files are loaded.
	ElectronPairProduction epp(CMB);
	epp.setPhotonField(IRB_Kneiske04);
	epp.setPhotonField(IRB_Stecker05);
	epp.setPhotonField(IRB_Franceschini08);
	epp.setPhotonField(IRB_Finke10);
	epp.setPhotonField(IRB_Dominguez11);
	epp.setPhotonField(IRB_Gilmore12);
	epp.setPhotonField(IRB_Stecker16_upper);
	epp.setPhotonField(IRB_Stecker16_lower);
}

TEST(ElectronPairProduction, energyDecreasing) {
	// Test if energy loss occurs for protons with energies from 1e15 - 1e23 eV.
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
	// Test if nothing happens below 1e15 eV.
	ElectronPairProduction epp(CMB);
	Candidate c(nucleusId(1, 1), 1E14 * eV);
	epp.process(&c);
	EXPECT_DOUBLE_EQ(1E14 * eV, c.current.getEnergy());
}

TEST(ElectronPairProduction, thisIsNotNucleonic) {
	// Test if non-nuclei are skipped.
	ElectronPairProduction epp(CMB);
	Candidate c(11, 1E20 * eV);  // electron
	epp.process(&c);
	EXPECT_DOUBLE_EQ(1E20 * eV, c.current.getEnergy());
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

// NuclearDecay ---------------------------------------------------------------
TEST(NuclearDecay, scandium44) {
	// Test beta+ decay of 44Sc to 44Ca.
	// This test can stochastically fail.
	NuclearDecay d(true, true);
	Candidate c(nucleusId(44, 21), 1E18 * eV);
	c.setCurrentStep(100 * Mpc);
	double gamma = c.current.getLorentzFactor();
	d.process(&c);
	
	// expected decay product: 44Ca
	EXPECT_EQ(nucleusId(44, 20), c.current.getId());

	// expect Lorentz factor to be conserved
	EXPECT_DOUBLE_EQ(gamma, c.current.getLorentzFactor());
	
	// expect at least two secondaries: positron + electron neutrino
	EXPECT_GE(c.secondaries.size(), 2);
}

TEST(NuclearDecay, lithium4) {
	// Test proton dripping of Li-4 to He-3
	// This test can stochastically fail
	NuclearDecay d;
	Candidate c(nucleusId(4, 3), 4 * EeV);
	c.setCurrentStep(100 * Mpc);
	d.process(&c);
	
	// expected decay product: He-3
	EXPECT_EQ(nucleusId(3, 2), c.current.getId());

	// expected secondary: proton
	EXPECT_EQ(1, c.secondaries.size());
	Candidate c1 = *c.secondaries[0];
	EXPECT_EQ(nucleusId(1, 1), c1.current.getId());
	EXPECT_EQ(1 * EeV, c1.current.getEnergy());
}

TEST(NuclearDecay, helium5) {
	// Test neutron dripping of He-5 to He-4.
	// This test can stochastically fail.
	NuclearDecay d;
	Candidate c(nucleusId(5, 2), 5 * EeV);
	c.setCurrentStep(100 * Mpc);
	d.process(&c);

	// expected primary: He-4
	EXPECT_EQ(nucleusId(4, 2), c.current.getId());
	EXPECT_EQ(4, c.current.getEnergy() / EeV);
	
	// expected secondary: neutron
	Candidate c2 = *c.secondaries[0];
	EXPECT_EQ(nucleusId(1, 0), c2.current.getId());
	EXPECT_EQ(1, c2.current.getEnergy() / EeV);
}

TEST(NuclearDecay, limitNextStep) {
	// Test if next step is limited in case of a neutron.
	NuclearDecay decay;
	Candidate c(nucleusId(1, 0), 10 * EeV);
	c.setNextStep(std::numeric_limits<double>::max());
	decay.process(&c);
	EXPECT_LT(c.getNextStep(), std::numeric_limits<double>::max());
}

TEST(NuclearDecay, allChannelsWorking) {
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

TEST(NuclearDecay, secondaries) {
	// Test if all types of secondaries are produced.
	NuclearDecay d;
	d.setHaveElectrons(true);
	d.setHaveNeutrinos(true);
	d.setHavePhotons(true);
	Candidate c;

	// He-8 --> Li-8 + e- + neutrino
	// additional photon emitted with 84% probability
	// --> expect at least 1 photon out of 10 decays
	for (int i = 0; i < 10; ++i) {
		c.current.setId(nucleusId(8, 2));
		c.current.setEnergy(5 * EeV);
		d.performInteraction(&c, 10000);
	}

	// count number of secondaries
	size_t nElectrons = 0;
	size_t nNeutrinos = 0;
	size_t nPhotons = 0;

	for(size_t i = 0; i < c.secondaries.size(); ++i) {
		int id = (*c.secondaries[i]).current.getId();
		if (id == 22) nPhotons++;
		if (id == 11) nElectrons++;
		if (id == -12) nNeutrinos++;
	}

	EXPECT_EQ(nElectrons, 10);
	EXPECT_EQ(nNeutrinos, 10);
	EXPECT_GE(nPhotons, 1);
}

TEST(NuclearDecay, thisIsNotNucleonic) {
	// Test if nothing happens to an electron
	NuclearDecay decay;
	Candidate c(11, 10 * EeV);
	c.setNextStep(std::numeric_limits<double>::max());
	decay.process(&c);
	EXPECT_EQ(11, c.current.getId());
	EXPECT_EQ(10 * EeV, c.current.getEnergy());
}

// PhotoDisintegration --------------------------------------------------------
TEST(PhotoDisintegration, allBackgrounds) {
	// Test if interaction data files are loaded.
	PhotoDisintegration pd(CMB);
	pd.setPhotonField(IRB_Kneiske04);
	pd.setPhotonField(IRB_Stecker05);
	pd.setPhotonField(IRB_Franceschini08);
	pd.setPhotonField(IRB_Finke10);
	pd.setPhotonField(IRB_Dominguez11);
	pd.setPhotonField(IRB_Gilmore12);
	pd.setPhotonField(IRB_Stecker16_upper);
	pd.setPhotonField(IRB_Stecker16_lower);
}

TEST(PhotoDisintegration, carbon) {
	// Test if a 100 EeV C-12 nucleus photo-disintegrates (at least once) over a distance of 1 Gpc.
	// This test can stochastically fail.
	PhotoDisintegration pd(CMB);
	Candidate c;
	int id = nucleusId(12, 6);
	c.current.setId(id);
	c.current.setEnergy(100 * EeV);
	c.setCurrentStep(1000 * Mpc);
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
	// Test if a 200 EeV Fe-56 nucleus photo-disintegrates (at least once) over a distance of 1 Gpc.
	// This test can stochastically fail.
	PhotoDisintegration pd(IRB);
	Candidate c;
	int id = nucleusId(56, 26);
	c.current.setId(id);
	c.current.setEnergy(200 * EeV);
	c.setCurrentStep(1000 * Mpc);
	pd.process(&c);

	// expect energy loss
	EXPECT_TRUE(c.current.getEnergy() < 200 * EeV);
	
	// expect secondaries produced
	EXPECT_TRUE(c.secondaries.size() > 0);

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

	// nucleon number conserved
	EXPECT_EQ(56, A);
	
	// proton number conserved (no decay active)
	EXPECT_EQ(26, Z);
	
	// energy conserved
	EXPECT_DOUBLE_EQ(200 * EeV, E);
}

TEST(PhotoDisintegration, thisIsNotNucleonic) {
	// Test that nothing happens to an electron.
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

TEST(Photodisintegration, updateParticleParentProperties) { // Issue: #204
	PhotoDisintegration pd(CMB);

	Candidate c(nucleusId(56,26), 500 * EeV, Vector3d(1 * Mpc, 0, 0));

	pd.performInteraction(&c, 1);
	// the candidates parent is the original particle
	EXPECT_EQ(c.created.getId(), nucleusId(56,26));

	pd.performInteraction(&c, 1);
	// now it has to be changed
	EXPECT_NE(c.created.getId(), nucleusId(56,26));
}


// ElasticScattering ----------------------------------------------------------
TEST(ElasticScattering, allBackgrounds) {
	// Test if interaction data files are loaded.
	ElasticScattering scattering(CMB);
	scattering.setPhotonField(IRB);
}

TEST(ElasticScattering, secondaries) {
	// Test the creation of cosmic ray photons.
	// This test can stochastically fail.
	ElasticScattering scattering(CMB);
	Candidate c;
	int id = nucleusId(12, 6);
	c.current.setId(id);
	c.current.setEnergy(200 * EeV);
	c.setCurrentStep(400 * Mpc);
	scattering.process(&c);

	EXPECT_GT(c.secondaries.size(), 0);

	for (int i = 0; i < c.secondaries.size(); i++) {
		int id = (*c.secondaries[i]).current.getId();
		EXPECT_EQ(id, 22);
		double energy = (*c.secondaries[i]).current.getEnergy();
		EXPECT_GT(energy, 0);
		EXPECT_LT(energy, 200 * EeV);
	}
}

// PhotoPionProduction --------------------------------------------------------
TEST(PhotoPionProduction, allBackgrounds) {
	// Test if all interaction data files can be loaded.
	PhotoPionProduction ppp;
	ppp.setPhotonField(IRB_Kneiske04);
	ppp.setPhotonField(IRB_Stecker05);
	ppp.setPhotonField(IRB_Franceschini08);
	ppp.setPhotonField(IRB_Finke10);
	ppp.setPhotonField(IRB_Dominguez11);
	ppp.setPhotonField(IRB_Gilmore12);
	ppp.setPhotonField(IRB_Stecker16_upper);
	ppp.setPhotonField(IRB_Stecker16_lower);
}

TEST(PhotoPionProduction, proton) {
	// Test photo-pion interaction for 100 EeV proton.
	// This test can stochastically fail.
	PhotoPionProduction ppp;
	Candidate c(nucleusId(1, 1), 100 * EeV);
	c.setCurrentStep(1000 * Mpc);
	ppp.process(&c);

	// expect energy loss
	EXPECT_LT(c.current.getEnergy(), 100 * EeV);

	// expect nucleon number conservation
	EXPECT_EQ(1, massNumber(c.current.getId()));

	// expect no (nucleonic) secondaries
	EXPECT_EQ(0, c.secondaries.size());
}

TEST(PhotoPionProduction, helium) {
	// Test photo-pion interaction for 400 EeV He nucleus.
	// This test can stochastically fail.
	PhotoPionProduction ppp;
	Candidate c;
	c.current.setId(nucleusId(4, 2));
	c.current.setEnergy(400 * EeV);
	c.setCurrentStep(1000 * Mpc);
	ppp.process(&c);
	EXPECT_LT(c.current.getEnergy(), 400 * EeV);
	int id = c.current.getId();
	EXPECT_TRUE(massNumber(id) < 4);
	EXPECT_TRUE(c.secondaries.size() > 0);
}

TEST(PhotoPionProduction, thisIsNotNucleonic) {
	// Test if noting happens to an electron.
	PhotoPionProduction ppp;
	Candidate c;
	c.current.setId(11); // electron
	c.current.setEnergy(10 * EeV);
	c.setCurrentStep(100 * Mpc);
	ppp.process(&c);
	EXPECT_EQ(11, c.current.getId());
	EXPECT_EQ(10 * EeV, c.current.getEnergy());
}

TEST(PhotoPionProduction, limitNextStep) {
	// Test if the interaction limits the next propagation step.
	PhotoPionProduction ppp;
	Candidate c(nucleusId(1, 1), 200 * EeV);
	c.setNextStep(std::numeric_limits<double>::max());
	ppp.process(&c);
	EXPECT_LT(c.getNextStep(), std::numeric_limits<double>::max());
}

TEST(PhotoPionProduction, secondaries) {
	// Test photo-pion interaction for 100 EeV proton.
	// This test can stochastically fail.
	PhotoPionProduction ppp(CMB, true, true, true);
	Candidate c(nucleusId(1, 1), 100 * EeV);
	c.setCurrentStep(1000 * Mpc);
	ppp.process(&c);
	// there should be secondaries
	EXPECT_GT(c.secondaries.size(), 1);
}

// Redshift -------------------------------------------------------------------
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

// EMPairProduction -----------------------------------------------------------
TEST(EMPairProduction, limitNextStep) {
	// Test if the interaction limits the next propagation step.
	EMPairProduction m;
	Candidate c(22, 1e17 * eV);
	c.setNextStep(std::numeric_limits<double>::max());
	m.process(&c);
	EXPECT_LT(c.getNextStep(), std::numeric_limits<double>::max());
}

TEST(EMPairProduction, secondaries) {
	// Test if secondaries are correctly produced.
	EMPairProduction m;
	m.setHaveElectrons(true);
	m.setThinning(0.);

	std::vector<PhotonField> fields;
	fields.push_back(CMB);
	fields.push_back(IRB_Gilmore12);
	fields.push_back(URB_Protheroe96);

	// loop over photon backgrounds
	for (int f = 0; f < fields.size(); f++) {
		m.setPhotonField(fields[f]);
		for (int i = 0; i < 130; i++) { // loop over energies Ep = (1e10 - 1e23) eV
			double Ep = pow(10, 10.05 + 0.1 * i) * eV;
			Candidate c(22, Ep);
			c.setCurrentStep(std::numeric_limits<double>::max());
			// c.setCurrentStep(1e10 * Mpc);
			m.process(&c);

			// pass if no interaction has occured (no tabulated rates)
			if (c.isActive())
				continue;
			
			// expect 2 secondaries
			EXPECT_EQ(c.secondaries.size(), 2);

			// expect electron / positron with energies 0 < E < Ephoton
			double Etot = 0;
			for (int j = 0; j < c.secondaries.size(); j++) {
				Candidate s = *c.secondaries[j];
				EXPECT_EQ(abs(s.current.getId()), 11);
				EXPECT_GT(s.current.getEnergy(), 0);
				EXPECT_LT(s.current.getEnergy(), Ep);
				Etot += s.current.getEnergy();
			}

			// test energy conservation
			EXPECT_DOUBLE_EQ(Ep, Etot);
		}
	}
}

// EMDoublePairProduction -----------------------------------------------------
TEST(EMDoublePairProduction, limitNextStep) {
	// Test if the interaction limits the next propagation step.
	EMDoublePairProduction m;
	Candidate c(22, 1E17 * eV);
	c.setNextStep(std::numeric_limits<double>::max());
	m.process(&c);
	EXPECT_LT(c.getNextStep(), std::numeric_limits<double>::max());
}

TEST(EMDoublePairProduction, secondaries) {
	// Test if secondaries are correctly produced.
	EMDoublePairProduction m;
	m.setHaveElectrons(true);
	m.setThinning(0.);

	std::vector<PhotonField> fields;
	fields.push_back(CMB);
	fields.push_back(IRB_Gilmore12);
	fields.push_back(URB_Protheroe96);

	// loop over photon backgrounds
	for (int f = 0; f < fields.size(); f++) {
		m.setPhotonField(fields[f]);
		
		// loop over energies Ep = (1E10 - 1E23) eV
		for (int i = 0; i < 130; i++) {
			double Ep = pow(10, 10.05 + 0.1 * i) * eV;
			Candidate c(22, Ep);
			c.setCurrentStep(std::numeric_limits<double>::max());
			m.process(&c);

			// pass if no interaction has occured (no tabulated rates)
			if (c.isActive())
				continue;
			
			// expect 2 secondaries (only one pair is considered)
			EXPECT_EQ(c.secondaries.size(), 2);

			// expect electron / positron with energies 0 < E < Ephoton
			double Etot = 0;
			for (int j = 0; j < c.secondaries.size(); j++) {
				Candidate s = *c.secondaries[j];
				EXPECT_EQ(abs(s.current.getId()), 11);
				EXPECT_GT(s.current.getEnergy(), 0);
				EXPECT_LT(s.current.getEnergy(), Ep);
				Etot += s.current.getEnergy();
			}

			// test energy conservation
			EXPECT_NEAR(Ep, Etot, 1E-9);
		}
	}
}

// EMTripletPairProduction ----------------------------------------------------
TEST(EMTripletPairProduction, limitNextStep) {
	// Test if the interaction limits the next propagation step.
	EMTripletPairProduction m;
	Candidate c(11, 1E17 * eV);
	c.setNextStep(std::numeric_limits<double>::max());
	m.process(&c);
	EXPECT_LT(c.getNextStep(), std::numeric_limits<double>::max());
}

TEST(EMTripletPairProduction, secondaries) {
	// Test if secondaries are correctly produced.
	EMTripletPairProduction m;
	m.setHaveElectrons(true);

	std::vector<PhotonField> fields;
	fields.push_back(CMB);
	fields.push_back(IRB_Gilmore12);
	fields.push_back(URB_Protheroe96);

	// loop over photon backgrounds
	for (int f = 0; f < fields.size(); f++) {
		m.setPhotonField(fields[f]);
		
		// loop over energies Ep = (1E10 - 1E23) eV
		for (int i = 0; i < 130; i++) {

			double Ep = pow(10, 10.05 + 0.1 * i) * eV;
			Candidate c(11, Ep);
			c.setCurrentStep(1e4 * Mpc); // use lower value so that the test can run faster
			m.process(&c);

			// pass if no interaction has occured (no tabulated rates)
			if (c.current.getEnergy() == Ep)
				continue;

			// expect positive energy of primary electron
			EXPECT_GT(c.current.getEnergy(), 0);
			double Etot = c.current.getEnergy();

			// expect electron / positron with energies 0 < E < Ephoton
			for (int j = 0; j < c.secondaries.size(); j++) {
				Candidate s = *c.secondaries[j];
				EXPECT_EQ(abs(s.current.getId()), 11);
				EXPECT_GT(s.current.getEnergy(), 0);
				EXPECT_LT(s.current.getEnergy(), Ep);
				Etot += s.current.getEnergy();
			}

			// test energy conservation
			EXPECT_NEAR(Ep, Etot, 1e-9);
		}
	}
}

// EMInverseComptonScattering -------------------------------------------------
TEST(EMInverseComptonScattering, limitNextStep) {
	// Test if the interaction limits the next propagation step.
	EMInverseComptonScattering m;
	Candidate c(11, 1E17 * eV);
	c.setNextStep(std::numeric_limits<double>::max());
	m.process(&c);
	EXPECT_LT(c.getNextStep(), std::numeric_limits<double>::max());
}

TEST(EMInverseComptonScattering, secondaries) {
	// Test if secondaries are correctly produced.
	EMInverseComptonScattering m;
	m.setHavePhotons(true);

	std::vector<PhotonField> fields;
	fields.push_back(CMB);
	fields.push_back(IRB_Gilmore12);
	fields.push_back(URB_Protheroe96);

	// loop over photon backgrounds
	for (int f = 0; f < fields.size(); f++) {
		m.setPhotonField(fields[f]);
		
		// loop over energies Ep = (1E10 - 1E23) eV
		for (int i = 0; i < 130; i++) {
			double Ep = pow(10, 10.05 + 0.1 * i) * eV;
			Candidate c(11, Ep);
			c.setCurrentStep(1e3 * Mpc); // use lower value so that the test can run faster
			m.process(&c);

			// pass if no interaction has occured (no tabulated rates)
			if (c.current.getEnergy() == Ep)
				continue;
			
			// expect positive energy of primary electron
			EXPECT_GT(c.current.getEnergy(), 0);

			// expect photon with energy 0 < E < Ephoton
			Candidate s = *c.secondaries[0];
			EXPECT_EQ(abs(s.current.getId()), 22);
			EXPECT_TRUE(s.current.getEnergy() >= 0.);
			EXPECT_TRUE(s.current.getEnergy() < Ep);


			double Etot = c.current.getEnergy();
			for (int j = 0; j < c.secondaries.size(); j++) {
				s = *c.secondaries[j];
				Etot += s.current.getEnergy();
			}
			EXPECT_NEAR(Ep, Etot, 1e-9); 
		}
	}
}


int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace crpropa

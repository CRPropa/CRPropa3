#include "crpropa/Units.h"
#include "crpropa/Common.h"
#include "crpropa/ParticleID.h"
#include "crpropa/module/CandidateSplitting.h"

#include "gtest/gtest.h"
#include <stdexcept>
#include <cmath>

namespace crpropa {

TEST(testCandidateSplitting, SimpleTest) {
	int nSplit = 2;
	int nBins = 4;
	double minWeight = pow(1. / nSplit, 2);
	double Emin = 1; // dimensionless for testing
	double Emax = 10; 	

	CandidateSplitting split_lin(nSplit, Emin, Emax, nBins, minWeight);
	double dE = (Emax - Emin) / nBins;
	EXPECT_DOUBLE_EQ(split_lin.getEnergyBins()[0], Emin);
	EXPECT_DOUBLE_EQ(split_lin.getEnergyBins()[1], Emin + dE);

	EXPECT_EQ(split_lin.getNsplit(), nSplit);
	EXPECT_DOUBLE_EQ(split_lin.getMinimalWeight(), minWeight);

	CandidateSplitting split_log(nSplit, Emin, Emax, nBins, minWeight, true);
	double dE_log = pow(Emax / Emin, 1. / (nBins - 1.0));
	EXPECT_DOUBLE_EQ(split_log.getEnergyBins()[0], Emin);
	EXPECT_DOUBLE_EQ(split_log.getEnergyBins()[1], Emin * dE_log);

	double spectralIndex = -2.;
	CandidateSplitting split_dsa(spectralIndex, Emin, nBins);
	double dE_dsa = pow(1. / 2, 1. / (spectralIndex + 1));
	EXPECT_DOUBLE_EQ(split_dsa.getEnergyBins()[0], Emin * dE_dsa);
	EXPECT_DOUBLE_EQ(split_dsa.getEnergyBins()[nBins - 1], Emin * pow(dE_dsa, nBins));
}


TEST(testCandidateSplitting, CheckSplits) {
	int nSplit = 2;
	int nBins = 3;
	double Emin = 1; // dimensionless for testing
	double Emax = 10;
	double minWeight = pow(1. / nSplit, 4);

	CandidateSplitting splitting(nSplit, Emin, Emax, nBins, minWeight);
	Candidate c(nucleusId(1,1),0.5);
	double weight = 1.0;
	double serial = c.getSerialNumber();
	
	splitting.process(&c); // no split
	EXPECT_DOUBLE_EQ(c.getWeight(), weight);
	EXPECT_DOUBLE_EQ(c.getNextSerialNumber(), serial);

	c.current.setEnergy(2); 
	splitting.process(&c); // 1. split
	weight = weight/nSplit;
	EXPECT_DOUBLE_EQ(c.getWeight(), weight);
	EXPECT_DOUBLE_EQ(c.getNextSerialNumber(), serial + 1);
	c.previous.setEnergy(2);

	c.current.setEnergy(6); 
	splitting.process(&c); // 2. split
	weight = weight/nSplit;
	EXPECT_DOUBLE_EQ(c.getWeight(), weight);
	EXPECT_DOUBLE_EQ(c.getNextSerialNumber(), serial + 2);
	c.previous.setEnergy(6);

	c.current.setEnergy(0.5); 
	splitting.process(&c); // no split, cooling
	EXPECT_DOUBLE_EQ(c.getWeight(), weight);
	EXPECT_DOUBLE_EQ(c.getNextSerialNumber(), serial + 2);
	c.previous.setEnergy(0.5);

	c.current.setEnergy(6); 
	splitting.process(&c); // 3. & 4. split, crosses two boundaries
	weight = weight/nSplit/nSplit;
	EXPECT_DOUBLE_EQ(c.getWeight(), weight);
	EXPECT_DOUBLE_EQ(c.getNextSerialNumber(), serial + 4);
	c.previous.setEnergy(6);

	c.current.setEnergy(8); 
	splitting.process(&c); // no split, minimal weight reached
	EXPECT_DOUBLE_EQ(c.getWeight(), weight);
	EXPECT_DOUBLE_EQ(c.getNextSerialNumber(), serial + 4);
	c.previous.setEnergy(8);
}

} //namespace crpropa

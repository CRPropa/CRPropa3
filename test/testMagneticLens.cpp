/// Tests gamale with a test lense from file and demonstrating basic
/// usage

//--------------------------------------------
// Project: Galactic magnetic Lense (GaMaLe) -
// Copyright (C) 2009 Tobias Winchen         -
//               RWTH Aachen, Germany        -
// Contact: winchen@physik.rwth-achen.de     -
// Licensed under the GNU GPL v2             - 
//--------------------------------------------

#include <string>
#include <boost/progress.hpp>
#include "gtest/gtest.h"

#include "parsec/MagneticLens.h"
#include "parsec/ModelMatrix.h"
#include "parsec/Pixelization.h"
#include "parsec/ParticleMapsContainer.h"

using namespace std;
using namespace parsec;

TEST(MagneticLens, Deflection)
{
	MagneticLens magneticLens(5);
	Pixelization P(5);
	ModelMatrix M(P.nPix(),P.nPix(),P.nPix());

	// Map any direction (p,t) to (p, -t)
	for (int i=0;i<P.nPix();i++)
	{
		double theta, phi;
		P.pix2Direction(i, phi, theta);
		theta*= -1;
		int j = P.direction2Pix(phi, theta);
		M(i,j) =1;
	}

	magneticLens.setLensPart(M,19.0,20.0);

	for (int i=0; i < P.nPix(); i++)
	{
		double theta, phi;
		P.pix2Direction(i, phi, theta);
		double theta0 = theta;

		// No CR is allowed to be lost in this lens
		EXPECT_TRUE(magneticLens.transformCosmicRay(19.32, phi, theta));
		EXPECT_NEAR(theta+theta0,0., 0.05);
	}
}

TEST(MagneticLens, OutOfBoundsEnergy)
{
	MagneticLens magneticLens(5);
	Pixelization P(5);
	ModelMatrix M(P.nPix(),P.nPix(),P.nPix());
	magneticLens.setLensPart(M,19.0,20.0);
	double theta, phi;
	EXPECT_THROW(magneticLens.transformCosmicRay(1.32, phi, theta), std::runtime_error);
}


TEST(Pixelization, angularDistance)
{
	// test for correct angular distance in case of same vectors 
	Pixelization P(6);
	for (int idx =0; idx < P.nPix(); idx++)
	{
		double ang = P.angularDistance(idx,idx);
		EXPECT_TRUE(ang == ang);
	}
}


TEST(ParticleMapsContainer, addParticle)
{
  ParticleMapsContainer maps;
  maps.addParticle(1000010010, 1E18, 0 , 0 );
  std::vector<int> pids = maps.getParticleIds();
  EXPECT_EQ(pids.size(), 1);
  EXPECT_EQ(pids[0], 1000010010);

  std::vector<double> energies = maps.getEnergies(1000010010);
  EXPECT_EQ(energies.size(), 1);
}

TEST(ParticleMapsContainer, getRandomParticles)
{
  ParticleMapsContainer maps(0.002);
  maps.addParticle(1000010010, 1E18, 0 , 0 );
  std::vector<double> energies;
  std::vector<double> lons;
  std::vector<double> lats;
  std::vector<int> particleIds;

  size_t N = 42;
	maps.getRandomParticles(N, particleIds, energies, lons, lats);
			
  EXPECT_EQ(energies.size(), N);
  EXPECT_EQ(lons.size(), N);
  EXPECT_EQ(lats.size(), N);
  EXPECT_EQ(particleIds.size(), N);

  for(size_t i = 0; i < N; i++)
  {
    EXPECT_NEAR(log10(energies[i]), 18, 0.002);
    EXPECT_EQ(particleIds[i], 1000010010);
    EXPECT_NEAR(lons[i], 0, 1./180*M_PI);
    EXPECT_NEAR(lats[i], 0, 1./180*M_PI);
  }

}

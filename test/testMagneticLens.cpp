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

#include "crpropa/magneticLens/MagneticLens.h"
#include "crpropa/magneticLens/ModelMatrix.h"
#include "crpropa/magneticLens/Pixelization.h"
#include "crpropa/magneticLens/ParticleMapsContainer.h"
#include "crpropa/Common.h"

using namespace std;
using namespace crpropa;

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

	magneticLens.setLensPart(M, 10 * EeV, 100 * EeV);

	for (int i=0; i < P.nPix(); i++)
	{
		double theta, phi;
		P.pix2Direction(i, phi, theta);
		double theta0 = theta;

		// No CR is allowed to be lost in this lens
		EXPECT_TRUE(magneticLens.transformCosmicRay(20 * EeV, phi, theta));
		EXPECT_NEAR(theta+theta0,0., 0.05);
	}
}

TEST(MagneticLens, Vector3Deflection)
{
	MagneticLens magneticLens(5);
	Pixelization P(5);
	ModelMatrix M(P.nPix(),P.nPix(),P.nPix());

	// No deflection 
	for (int i=0;i<P.nPix();i++)
	{
		M(i,i) = 1;
	}

	magneticLens.setLensPart(M, 10 * EeV, 100 * EeV);

	Vector3d u(1,0,0);
	Vector3d v(u);
	EXPECT_TRUE(magneticLens.transformCosmicRay(20 * EeV, v));
	EXPECT_NEAR(u.getAngleTo(v), 0., 2. / 180 * M_PI);

	u.x = 0; u.y = 1; u.z = 0; v = u;
	EXPECT_TRUE(magneticLens.transformCosmicRay(20 * EeV, v));
	EXPECT_NEAR(u.getAngleTo(v), 0., 2. / 180 * M_PI);

	u.x = 0; u.y = 0; u.z = 1; v = u;
	EXPECT_TRUE(magneticLens.transformCosmicRay(20 * EeV, v));
	EXPECT_NEAR(u.getAngleTo(v), 0., 2. / 180 * M_PI);

	u.x = 1; u.y = -2; u.z = 3; v = u;
	EXPECT_TRUE(magneticLens.transformCosmicRay(20 * EeV, v));
	EXPECT_NEAR(u.getAngleTo(v), 0., 2. / 180 * M_PI);
}

TEST(MagneticLens, OutOfBoundsEnergy)
{
	MagneticLens magneticLens(5);
	Pixelization P(5);
	ModelMatrix M(P.nPix(),P.nPix(),P.nPix());
	magneticLens.setLensPart(M,10. * EeV, 100. * EeV);
	double theta, phi;
	EXPECT_FALSE(magneticLens.transformCosmicRay(1. * EeV, phi, theta));
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
  maps.addParticle(1000010010, 1 * EeV, 0 , 0 );
  std::vector<int> pids = maps.getParticleIds();
  EXPECT_EQ(pids.size(), 1);
  EXPECT_EQ(pids[0], 1000010010);

  std::vector<double> energies = maps.getEnergies(1000010010);
  EXPECT_EQ(energies.size(), 1);
}

TEST(ParticleMapsContainer, getRandomParticles)
{
  ParticleMapsContainer maps(0.002);
  maps.addParticle(1000010010, 1 * EeV, 0 , 0 );
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

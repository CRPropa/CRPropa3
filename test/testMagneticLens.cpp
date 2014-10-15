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

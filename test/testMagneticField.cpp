#include <stdexcept>

#include "crpropa/magneticField/MagneticFieldGrid.h"
#include "crpropa/magneticField/CMZField.h"
#include "crpropa/magneticField/PolarizedSingleModeMagneticField.h"
#include "crpropa/magneticField/GalacticMagneticField.h"
#include "crpropa/magneticField/UF23Field.h"
#include "crpropa/Grid.h"
#include "crpropa/Units.h"
#include "crpropa/Common.h"

#include "gtest/gtest.h"

using namespace crpropa;

TEST(testUniformMagneticField, SimpleTest) {
	UniformMagneticField B(Vector3d(-1, 5, 3));
	Vector3d b = B.getField(Vector3d(1, 0, 0));
	EXPECT_DOUBLE_EQ(b.getX(), -1);
	EXPECT_DOUBLE_EQ(b.getY(), 5);
	EXPECT_DOUBLE_EQ(b.getZ(), 3);
}

TEST(testMagneticDipoleField, SimpleTest) {
	// Test magnetic dipole
	// mu0 / (4*M_PI) * m / r^3 (2*cos(theta)*e_r + sin(theta)*e_theta)
	MagneticDipoleField B(Vector3d(0,0,0), Vector3d(0,0,1), 1);
	Vector3d b1 = B.getField(Vector3d(0, 0, 1)); // theta = 0
	Vector3d b2 = B.getField(Vector3d(1, 0, 0)); // theta = 0
	EXPECT_NEAR(b1.getX(), 0, 1E-8);
	EXPECT_NEAR(b1.getY(), 0, 1E-8);
	EXPECT_NEAR(b1.getZ(), mu0 / (4*M_PI) * 2, 1E-8);
	EXPECT_NEAR(b2.getX(), 0, 1E-8);
	EXPECT_NEAR(b2.getY(), 0, 1E-8);
	EXPECT_NEAR(b2.getZ(), -1 * mu0 / (4*M_PI), 1E-8);
}

#ifdef CRPROPA_HAVE_MUPARSER
TEST(testRenormalizeMagneticField, simpleTest) {
	ref_ptr<UniformMagneticField> field = new UniformMagneticField(Vector3d(2*nG, 0, 0));
	RenormalizeMagneticField modField(field, "B^2-1*nG");
	Vector3d b = modField.getField(Vector3d(5, 5, 5));
	EXPECT_NEAR(b.getR(), 3*nG, 0.001);
}
#endif

TEST(testMagneticFieldList, SimpleTest) {
	// Test a list of three magnetic fields
	MagneticFieldList B;
	B.addField(new UniformMagneticField(Vector3d(1, 0, 0)));
	B.addField(new UniformMagneticField(Vector3d(0, 2, 0)));
	B.addField(new UniformMagneticField(Vector3d(0, 0, 3)));
	Vector3d b = B.getField(Vector3d(0.));
	EXPECT_DOUBLE_EQ(b.x, 1);
	EXPECT_DOUBLE_EQ(b.y, 2);
	EXPECT_DOUBLE_EQ(b.z, 3);
}

TEST(testMagneticFieldEvolution, SimpleTest) {
	// Test if this decorator scales the underlying field as (1+z)^m
	ref_ptr<UniformMagneticField> B = new UniformMagneticField(Vector3d(1,0,0));
	double z = 1.2;
	double m = 3;
	MagneticFieldEvolution Bz(B, m);

	// scaled field
	Vector3d b = Bz.getField(Vector3d(0,0,0), z);
	EXPECT_DOUBLE_EQ(b.x, pow(1+z, m));

	// unscaled field
	b = Bz.getField(Vector3d(0,0,0));
	EXPECT_DOUBLE_EQ(b.x, 1);
}

class EchoMagneticField: public MagneticField {
public:
	Vector3d getField(const Vector3d &position) const {
		return position;
	}
};

TEST(testPeriodicMagneticField, Exceptions) {
	ref_ptr<EchoMagneticField> f = new EchoMagneticField();
	ref_ptr<PeriodicMagneticField> p = new PeriodicMagneticField(f,
			Vector3d(10000, 10000, 10000), Vector3d(1000, 1000, 1000), true);

	// box 0, 0, 0
	Vector3d v = p->getField(Vector3d(1000, 2000, 3000));
	EXPECT_DOUBLE_EQ(0, v.x);
	EXPECT_DOUBLE_EQ(1000, v.y);
	EXPECT_DOUBLE_EQ(2000, v.z);

	// box 1, 2, 3
	v = p->getField(Vector3d(12000, 23000, 34000));
	EXPECT_DOUBLE_EQ(9000, v.x);
	EXPECT_DOUBLE_EQ(2000, v.y);
	EXPECT_DOUBLE_EQ(7000, v.z);

	// box -1, -2, -3
	v = p->getField(Vector3d(0, -10000, -20000));
	EXPECT_DOUBLE_EQ(1000, v.x);
	EXPECT_DOUBLE_EQ(9000, v.y);
	EXPECT_DOUBLE_EQ(1000, v.z);

}

TEST(testCMZMagneticField, SimpleTest) {
	ref_ptr<CMZField> field = new CMZField();

	// check use-Values
	EXPECT_FALSE(field->getUseMCField());
	EXPECT_TRUE(field->getUseICField());
	EXPECT_FALSE(field->getUseNTFField());
	EXPECT_FALSE(field->getUseRadioArc());

	// check set function
	field->setUseMCField(true);
	EXPECT_TRUE(field->getUseMCField());
	field->setUseICField(false);
	EXPECT_FALSE(field->getUseICField());
	field->setUseNTFField(true);
	EXPECT_TRUE(field->getUseNTFField());
	field->setUseRadioArc(true);
	EXPECT_TRUE(field->getUseRadioArc());
}

TEST(testCMZMagneticField, TestICComponent) {
	ref_ptr<CMZField> field = new CMZField();
	Vector3d pos(10*pc,15*pc,-5*pc);

	// check IC field at given position
	Vector3d bVec = field->getField(pos);
	EXPECT_NEAR(bVec.getR(), 10.501*muG, 1E-3*muG);
	EXPECT_NEAR(bVec.x, 0.225*muG, 1E-3*muG);
	EXPECT_NEAR(bVec.y, 0.524*muG, 1E-3*muG);
	EXPECT_NEAR(bVec.z, 10.486*muG, 1E-3*muG);
}
TEST(testCMZMagneticField, TestNTFField){
	ref_ptr<CMZField> field = new CMZField();
	Vector3d pos(10*pc,15*pc,-5*pc);

	// check NFTField at given position
	Vector3d bVec = field->getNTFField(pos);
	EXPECT_NEAR(bVec.getR(),1.692*muG, 1e-3*muG);
	EXPECT_NEAR(bVec.x, -0.584*muG, 1e-3*muG);
	EXPECT_NEAR(bVec.y, -1.185*muG, 1e-3*muG);
	EXPECT_NEAR(bVec.z, 1.057*muG, 1e-3*muG);
}
TEST(testCMZMagneticField, TestRaidoArcField){
	ref_ptr<CMZField> field = new CMZField();
	Vector3d pos(10*pc,15*pc,-5*pc);

	// check RadioArcField at given position
	Vector3d bVec = field->getRadioArcField(pos);
	EXPECT_NEAR(bVec.getR(), 31.616*muG, 1e-3*muG);
	EXPECT_NEAR(bVec.x, -4.671*muG, 1e-3*muG);
	EXPECT_NEAR(bVec.y, 5.465*muG, 1e-3*muG);
	EXPECT_NEAR(bVec.z, 30.788*muG, 1e-3*muG);
}

TEST(testCMZMagneticField, TestAzimutalComponent){
	ref_ptr<CMZField> field = new CMZField();
	Vector3d mid(12*pc, 9*pc, 20*pc);
	Vector3d pos(9*pc, 10*pc, 25*pc);

	// simple Test for inner part
	Vector3d bVec = field->BAz(pos, mid, 100, 0.2, 60*pc);
	EXPECT_NEAR(bVec.x, 3939.782, 1e-3);
	EXPECT_NEAR(bVec.y, 14347.304, 1e-3);
	EXPECT_DOUBLE_EQ(bVec.z, 0);

	// simple Test for outer part
	bVec = field->BAz(pos, mid, 100, 0.2, 10*pc);
	EXPECT_NEAR(bVec.x, -164.659, 1e-3);
	EXPECT_NEAR(bVec.y, -1317.270, 1e-3);
	EXPECT_DOUBLE_EQ(bVec.z, 0);

	// test for molecular Clouds
	bVec = field->getMCField(pos);
	EXPECT_NEAR(bVec.x, -8.339*muG, 1e-3*muG);
	EXPECT_NEAR(bVec.y, -0.850*muG, 1e-3*muG);
	EXPECT_DOUBLE_EQ(bVec.z, 0);
}

TEST(testToroidalHaloField, SimpleTest) {
	ref_ptr<ToroidalHaloField> field = new ToroidalHaloField();
	Vector3d b = field->getField(Vector3d(0.));
	EXPECT_DOUBLE_EQ(b.x, 0);
	EXPECT_DOUBLE_EQ(b.y, 0);
	EXPECT_DOUBLE_EQ(b.z, 0);

	b = field->getField(Vector3d(1,0,0));
	EXPECT_DOUBLE_EQ(b.x, 0.5);
	EXPECT_DOUBLE_EQ(b.y, 0);
	EXPECT_DOUBLE_EQ(b.z, 0);
}

TEST(testLogarithmicSpiralField, SimpleTest) {
	ref_ptr<LogarithmicSpiralField> field = new LogarithmicSpiralField();
	Vector3d b = field->getField(Vector3d(8.5, 0, 0)*kpc);
	EXPECT_NEAR(b.x, -1., 1E-2);
	EXPECT_NEAR(b.y, 0, 1E-10);
	EXPECT_NEAR(b.z, 0, 1E-10);
}


// UF23 reference values
std::vector<std::vector<Vector3d>>
getUF23ReferenceValues()
{
  using V = Vector3d;
  std::vector<std::vector<Vector3d>> retVal;
  retVal.push_back({
                    V{-1.46860417e+00, 1.64555489e+00, 8.35702311e-01},
                    V{0.00000000e+00, -4.22433497e-01, 1.83232985e-01},
                    V{-3.00239177e-01, -2.97045767e-01, 1.83232985e-01},
                    V{8.58382023e-04, 7.71891049e-03, 9.71705527e-01},
                    V{-1.17276875e+00, -2.33013590e-01, 4.10798494e-01},
                    V{2.63883569e-01, -8.79631081e-01, 5.54229961e-06},
                    V{-1.71971907e-02, -2.00291358e-02, 4.16024875e-02},
                    V{6.18960813e-01, -2.14437026e+00, 8.35702315e-01},
                    V{-1.50883911e+00, 2.63874909e+00, 1.16578928e-02}
    });
  retVal.push_back({
                    V{-1.39447625e+00, 1.57135437e+00, 8.65163870e-01},
                    V{0.00000000e+00, -4.44318027e-01, 1.54430718e-01},
                    V{-3.15695875e-01, -3.12652066e-01, 1.54430718e-01},
                    V{2.09678202e-03, 2.63832916e-03, 9.81077691e-01},
                    V{-1.08603376e+00, -1.86359136e-01, 3.91730436e-01},
                    V{2.11095801e-01, -7.04082667e-01, 8.86591950e-05},
                    V{-1.38087843e-02, -1.52011002e-02, 3.06842086e-02},
                    V{5.76543516e-01, -2.03544911e+00, 8.65163870e-01},
                    V{-3.77778745e-01, 6.89644787e-01, 5.85347897e-02},
    });
  retVal.push_back({
                    V{-1.04170166e+00, 1.93212739e+00, 2.95362499e+00},
                    V{0.00000000e+00, -5.87995638e-01, 4.66925506e-01},
                    V{-4.20802701e-01, -4.10291633e-01, 4.59557272e-01},
                    V{7.79982996e-03, 1.50482351e-02, 5.50337531e+00},
                    V{-1.34639246e+00, 6.80103301e-02, 8.60002649e-01},
                    V{1.37242004e-01, -8.67085225e-01, 1.25106474e-01},
                    V{-1.69325979e-02, -2.95454738e-02, 4.72616404e-02},
                    V{4.53873303e-01, -2.03292501e+00, 9.95205454e-01},
                    V{-1.30599234e+00, 2.25151057e+00, 2.56704152e-01},
    });
  retVal.push_back({
                    V{-1.45840221e+00, 1.64227817e+00, 8.41586777e-01},
                    V{0.00000000e+00, -6.98341816e-01, 1.84323895e-01},
                    V{-4.95342500e-01, -4.92154113e-01, 1.84323895e-01},
                    V{-5.27533653e-03, 1.48965596e-02, 9.86107885e-01},
                    V{-1.47472798e+00, -3.33905274e-01, 4.11265918e-01},
                    V{2.31202234e-01, -7.70722755e-01, 1.12129423e-05},
                    V{-1.54200071e-02, -2.22012150e-02, 4.22733003e-02},
                    V{5.92201401e-01, -2.10378256e+00, 8.41586782e-01},
                    V{4.05057122e-03, -9.76673905e-03, 8.24321504e-03},
    });
  retVal.push_back({
                    V{-1.50595160e+00, 1.67869437e+00, 8.29806168e-01},
                    V{0.00000000e+00, -2.27289811e-01, 1.86869632e-01},
                    V{-1.62291024e-01, -1.59076795e-01, 1.86869632e-01},
                    V{6.71536463e-04, 7.88990042e-03, 9.63291786e-01},
                    V{-9.35561196e-01, -1.55840707e-01, 4.13609050e-01},
                    V{2.68120686e-01, -8.93743434e-01, 3.48927149e-06},
                    V{-1.88291055e-02, -1.94245715e-02, 4.29970299e-02},
                    V{6.50804557e-01, -2.18817590e+00, 8.29806172e-01},
                    V{-1.58450595e+00, 2.77817014e+00, 1.10336683e-02},
    });
  retVal.push_back({
                    V{-1.33098823e+00, 1.49556466e+00, 6.88642986e-01},
                    V{0.00000000e+00, -4.99342453e-01, 1.16143813e-01},
                    V{-3.54225500e-01, -3.51947350e-01, 1.16143813e-01},
                    V{2.14352986e-03, 3.57971370e-03, 8.05219599e-01},
                    V{-1.15136627e+00, -2.48536546e-01, 2.95713475e-01},
                    V{1.18057505e-01, -3.98308301e-01, 9.11299352e-04},
                    V{-1.06899006e-02, -1.12335952e-02, 2.33358311e-02},
                    V{5.61264936e-01, -1.94601963e+00, 6.88642987e-01},
                    V{-1.18200258e+00, 2.01179003e+00, 8.02453884e-02},
    });
  retVal.push_back({
                    V{-3.65508814e-01, 4.74456797e-01, 5.76165841e-01},
                    V{0.00000000e+00, 0.00000000e+00, 6.36668108e-02},
                    V{-3.68784336e-03, -2.20689056e-03, 6.36668108e-02},
                    V{-5.03208784e-02, 5.09887480e-02, 6.27665021e-01},
                    V{-4.04821314e-01, -1.01426396e-02, 2.06022291e-01},
                    V{-2.34287239e-02, -9.72536117e-02, 2.80624948e-02},
                    V{-4.42698169e-03, -6.23429869e-03, 1.07555956e-02},
                    V{3.25022555e-01, -1.18950549e+00, 5.76163734e-01},
                    V{-8.65217092e-01, 1.85270104e+00, 3.73913202e-01},
    });
  retVal.push_back({
                    V{-1.92580231e+00, 2.17134243e+00, 1.13723929e+00},
                    V{0.00000000e+00, -4.73261099e-01, 2.65974345e-01},
                    V{-3.36780079e-01, -3.32368900e-01, 2.65974345e-01},
                    V{-5.24558924e-04, 1.51886238e-02, 1.33757258e+00},
                    V{-1.47213976e+00, -2.82203258e-01, 5.71920142e-01},
                    V{3.54206333e-01, -1.18108962e+00, 1.00077423e-04},
                    V{-2.67268842e-02, -2.88588605e-02, 6.38364561e-02},
                    V{7.90570955e-01, -2.84039611e+00, 1.13723931e+00},
                    V{-1.97316856e+00, 3.44153600e+00, 3.20266742e-02},
    });
  return retVal;
}


TEST(testUF23Field, SimpleTest) {
  const std::vector<UF23Field::ModelType> models =
    {
     UF23Field::base, UF23Field::neCL, UF23Field::expX, UF23Field::spur,
     UF23Field::cre10, UF23Field::synCG, UF23Field::twistX, UF23Field::nebCor
    };

  using V = Vector3d;
  const std::vector<Vector3d> testPositions =
    {
     V{1, 1, 1}, V{0, 0, -8}, V{0.1, -0.1, -8}, V{0.1, 0.1, 0.1}, V{-1, 3, 4},
     V{-10, -3, 2}, V{-10, -10, 20}, V{-4, -2, 1}, V{6, 5, -0.1}
    };

  const std::vector<std::vector<Vector3d>> referenceValues =
    getUF23ReferenceValues();

  const double precision = 1e-6;

  for (unsigned int i = 0; i < models.size(); ++i) {
    const auto& model = models[i];
    const UF23Field uf23Field(model);
    for (unsigned int j = 0; j < testPositions.size(); ++j) {
      const auto& position = testPositions[j];
      const auto val = uf23Field.getField(position*kpc);
      const auto& refVal = referenceValues[i][j] * microgauss;
      EXPECT_NEAR(val.x, refVal.x, precision);
      EXPECT_NEAR(val.y, refVal.y, precision);
      EXPECT_NEAR(val.z, refVal.z, precision);
    }
  }
}

TEST(testPolarizedSingleModeMagneticField, SimpleTest) {
	PolarizedSingleModeMagneticField B(2, 4, 0.5, Vector3d(1,1,1), Vector3d(0,1,0), Vector3d(1,0,0), "amplitude", "polarization", "elliptical");
	Vector3d b = B.getField(Vector3d(1,1,2));
	EXPECT_DOUBLE_EQ(b.x, 1);
	EXPECT_NEAR(b.y, 0, 1E-10);
	EXPECT_NEAR(b.z, 0, 1E-10);
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

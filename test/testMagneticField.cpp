#include <stdexcept>

#include "crpropa/magneticField/MagneticFieldGrid.h"
#include "crpropa/magneticField/CMZField.h"
#include "crpropa/magneticField/PolarizedSingleModeMagneticField.h"
#include "crpropa/magneticField/GalacticMagneticField.h"
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

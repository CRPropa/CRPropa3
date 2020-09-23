#include <stdexcept>

#include "crpropa/magneticField/MagneticFieldGrid.h"
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

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

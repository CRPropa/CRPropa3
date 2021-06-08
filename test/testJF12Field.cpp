/** Unit tests for Output modules of CRPropa
    JF12Field
 */

#include "CRPropa.h"
#include "gtest/gtest.h"

namespace crpropa {

// These tests just check consistency across changes in the JF12 code
// implemented in CRPropa3,
// one should compare these with the values from the JF12 paper to check
// the consistency with the described model
TEST(JF12Field, regularComponent) {
    auto JF12 = JF12Field();
    auto B1 = JF12.getRegularField(Vector3d(1*kpc, 0, 0));
    EXPECT_DOUBLE_EQ(B1.getX(), 0);
    EXPECT_DOUBLE_EQ(B1.getY(), 0);
    EXPECT_DOUBLE_EQ(B1.getZ(), 0);
    auto B2 = JF12.getRegularField(Vector3d(8*kpc, 0, 0));
    EXPECT_NEAR(B2.getR(), 1.7*muG, 0.1*muG);
}

TEST(JF12Field, turbulentComponent) {
    auto JF12 = JF12Field();
    JF12.randomTurbulent(1);
    auto B1 = JF12.getTurbulentField(Vector3d(1*kpc, 0, 0));
    EXPECT_NEAR(B1.getR(), 10.9*muG, 0.1*muG);
    auto B2 = JF12.getTurbulentField(Vector3d(8*kpc, 0, 0));
    EXPECT_NEAR(B2.getR(), 3.64*muG, 0.01*muG);
    auto B3 = JF12.getTurbulentField(Vector3d(3*kpc, 3*kpc, 0.1*kpc));
    EXPECT_NEAR(B3.getR(), 11.5*muG, 0.1*muG);
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace crpropa

#include "mpc/module/common.h"

#include "gtest/gtest.h"

namespace mpc {

TEST(common, digit) {
	EXPECT_EQ(1, digit(1234, 1000));
	EXPECT_EQ(2, digit(1234, 100));
	EXPECT_EQ(3, digit(1234, 10));
	EXPECT_EQ(4, digit(1234, 1));
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace mpc

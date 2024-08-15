#include "P528GTest.h"

class TestP528 : public ::testing::Test {
protected:
    void SetUp() override {
        // Load test data from CSV
        testData = ReadP528InputsAndResultTestData("P528Test_001.csv");
    }

    // Vector to hold test data
    std::vector<InputsAndResult> testData;

    double TOLERANCE = 1.0e-10;
};

// Test case to verify P528 Success
TEST_F(TestP528, TestP528Success) {
    // Ensure test data was loaded
    EXPECT_NE(static_cast<int>(testData.size()), 0);

    for (const auto& data : testData) {
        Result result;
        int rtn = P528(data.d__km, data.h_1__meter, data.h_2__meter, data.f__mhz, data.T_pol, data.p, &result);

        EXPECT_EQ(rtn, data.expectedReturn);
        EXPECT_EQ(result.propagation_mode, data.expectedResult.propagation_mode);
        EXPECT_EQ(result.warnings, data.expectedResult.warnings);
        EXPECT_NEAR(result.d__km, data.expectedResult.d__km, 1.0e-3);
        EXPECT_NEAR(result.A__db, data.expectedResult.A__db, TOLERANCE);
        EXPECT_NEAR(result.A_fs__db, data.expectedResult.A_fs__db, TOLERANCE);
        EXPECT_NEAR(result.A_a__db, data.expectedResult.A_a__db, TOLERANCE);
        EXPECT_NEAR(result.theta_h1__rad, data.expectedResult.theta_h1__rad, TOLERANCE);
    }
}

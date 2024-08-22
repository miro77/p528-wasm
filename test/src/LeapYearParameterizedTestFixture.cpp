#include "P528GTest.h"

class LeapYearParameterizedTestFixture :public ::testing::TestWithParam<InputsAndResult> {
protected:
    
};

TEST_P(LeapYearParameterizedTestFixture, OddYearsAreNotLeapYears) {
    InputsAndResult data = GetParam();
    ASSERT_EQ(data.d__km, 0);
}

INSTANTIATE_TEST_CASE_P(
    LeapYearTests,
    LeapYearParameterizedTestFixture,
    ::testing::ValuesIn(getAllFilesInDirectory())
);



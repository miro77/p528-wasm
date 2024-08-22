#include "P528GTest.h"
#include <cctype>
#include <string>
#include <functional>

class TestP528DataTables : public ::testing::TestWithParam<std::string>
{
public:
    struct PrintToStringParamName
    {
        template <class ParamType>
        std::string operator() (const testing::TestParamInfo<ParamType>& info) const
        {
            std::string testName = static_cast<std::string>(info.param);
            testName.resize(testName.find_last_of("."));
            testName.erase(remove_if(testName.begin(), testName.end(), [](char c) { return !isalnum(c); }), testName.end());

            return testName;
        }
    };

    double TOLERANCE = 1.0e-1;
    int testStep = 100;
};

TEST_P(TestP528DataTables, TestDataTable) {
    std::string fileName = GetParam();
    std::vector<InputsAndResult> testData = ReadP528InputsAndResultFromDataTable(fileName, testStep);
    EXPECT_NE(static_cast<int>(testData.size()), 0);
    std::cout << "TestP528 from Data Table '" << fileName << "', " << testData.size() << " Test instances." << std::endl;
    for (const auto& data : testData) {
        Result result;
        int rtn = P528(data.d__km, data.h_1__meter, data.h_2__meter, data.f__mhz, data.T_pol, data.p, &result);
        if (data.expectedResult.A_fs__db != NULL) {
            EXPECT_NEAR(result.A_fs__db, data.expectedResult.A_fs__db, TOLERANCE);
        }
        EXPECT_NEAR(result.A__db, data.expectedResult.A__db, TOLERANCE);
    }
}

INSTANTIATE_TEST_CASE_P(
    TestP528, 
    TestP528DataTables, 
    ::testing::ValuesIn(GetDataTablesFileList()));
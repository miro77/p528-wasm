#include "P528GTest.h"
#include "p528.h"

#include <vector>

using namespace ITS::ITU::PSeries::p528;

/*=============================================================================
 |
 |  Description:  Tests for P528_Profile(). The point-by-point function
 |                P528() is the reference: a profile must return the same
 |                result for every distance. These tests need no data files.
 |
 *===========================================================================*/
class TestP528Profile : public ::testing::Test {
protected:
    // The profile interpolates the gaseous attenuation, which changes the
    // results by rounding-level amounts only.
    double TOLERANCE = 1.0e-6;

    // distances from 0 km to far beyond the horizon
    std::vector<double> MakeDistances(double d_max__km, int n) {
        std::vector<double> d(n);
        for (int i = 0; i < n; i++)
            d[i] = d_max__km * i / (n - 1);
        return d;
    }

    void ExpectProfileMatchesPointwise(const std::vector<double>& d, double h_1__meter,
        double h_2__meter, double f__mhz, int T_pol, double p) {
        std::vector<Result> results(d.size());
        std::vector<int> return_codes(d.size());
        P528_Profile(d.data(), static_cast<int>(d.size()), h_1__meter, h_2__meter, f__mhz,
            T_pol, p, results.data(), return_codes.data());

        for (size_t i = 0; i < d.size(); i++) {
            SCOPED_TRACE("d__km = " + std::to_string(d[i]));
            Result expected;
            int rtn = P528(d[i], h_1__meter, h_2__meter, f__mhz, T_pol, p, &expected);

            EXPECT_EQ(return_codes[i], rtn);
            EXPECT_EQ(results[i].propagation_mode, expected.propagation_mode);
            EXPECT_EQ(results[i].warnings, expected.warnings);
            EXPECT_NEAR(results[i].d__km, expected.d__km, TOLERANCE);
            EXPECT_NEAR(results[i].A__db, expected.A__db, TOLERANCE);
            EXPECT_NEAR(results[i].A_fs__db, expected.A_fs__db, TOLERANCE);
            EXPECT_NEAR(results[i].A_a__db, expected.A_a__db, TOLERANCE);
            EXPECT_NEAR(results[i].theta_h1__rad, expected.theta_h1__rad, TOLERANCE);
        }
    }
};

TEST_F(TestP528Profile, MatchesPointwiseGroundToAir) {
    ExpectProfileMatchesPointwise(MakeDistances(600, 25), 10, 1000, 1200,
        POLARIZATION__HORIZONTAL, 50);
}

TEST_F(TestP528Profile, MatchesPointwiseAirToAir) {
    // an elevated low terminal gives downward-looking rays near the horizon
    ExpectProfileMatchesPointwise(MakeDistances(1200, 25), 1000, 12000, 5100,
        POLARIZATION__VERTICAL, 5);
}

TEST_F(TestP528Profile, MatchesPointwiseLowTerminals) {
    ExpectProfileMatchesPointwise(MakeDistances(300, 25), 5, 30, 450,
        POLARIZATION__VERTICAL, 95);
}

TEST_F(TestP528Profile, MatchesPointwiseWithHeightWarning) {
    // terminals above 20 km raise warnings for every point
    ExpectProfileMatchesPointwise(MakeDistances(1500, 15), 3000, 25000, 15500,
        POLARIZATION__HORIZONTAL, 50);
}

TEST_F(TestP528Profile, DistancesInAnyOrder) {
    // beyond the horizon first, then line-of-sight, then beyond again
    ExpectProfileMatchesPointwise({ 400, 20, 0, 250, 130.5, 20, 800 }, 10, 1000, 1200,
        POLARIZATION__HORIZONTAL, 50);
}

TEST_F(TestP528Profile, SameHeightsAtZeroDistance) {
    // the terminals coincide for d = 0, P528() reports that as SUCCESS
    ExpectProfileMatchesPointwise({ 0, 10, 0, 300 }, 100, 100, 1200,
        POLARIZATION__HORIZONTAL, 50);
}

TEST_F(TestP528Profile, InvalidDistanceOnlyAffectsItsPoint) {
    std::vector<double> d = { 50, -1, 400 };
    std::vector<Result> results(d.size());
    std::vector<int> return_codes(d.size());
    int rtn = P528_Profile(d.data(), 3, 10, 1000, 1200, POLARIZATION__HORIZONTAL, 50,
        results.data(), return_codes.data());

    EXPECT_EQ(rtn, ERROR_VALIDATION__D_KM);
    EXPECT_EQ(return_codes[0], SUCCESS);
    EXPECT_EQ(return_codes[1], ERROR_VALIDATION__D_KM);
    EXPECT_EQ(return_codes[2], SUCCESS);
    EXPECT_EQ(results[1].propagation_mode, PROP_MODE__NOT_SET);
    EXPECT_GT(results[0].A__db, 0);
    EXPECT_GT(results[2].A__db, results[0].A__db);
}

TEST_F(TestP528Profile, InvalidInputsFailEveryPoint) {
    std::vector<double> d = { 50, 100 };
    std::vector<Result> results(d.size());
    std::vector<int> return_codes(d.size());

    EXPECT_EQ(P528_Profile(d.data(), 2, 1, 1000, 1200, POLARIZATION__HORIZONTAL, 50,
        results.data(), return_codes.data()), ERROR_VALIDATION__H_1);
    EXPECT_EQ(return_codes[0], ERROR_VALIDATION__H_1);
    EXPECT_EQ(return_codes[1], ERROR_VALIDATION__H_1);

    EXPECT_EQ(P528_Profile(d.data(), 2, 10, 1000, 50, POLARIZATION__HORIZONTAL, 50,
        results.data(), return_codes.data()), ERROR_VALIDATION__F_MHZ_LOW);
    EXPECT_EQ(P528_Profile(d.data(), 2, 10, 1000, 1200, 2, 50,
        results.data(), return_codes.data()), ERROR_VALIDATION__POLARIZATION);
    EXPECT_EQ(P528_Profile(d.data(), 2, 10, 1000, 1200, POLARIZATION__HORIZONTAL, 0.5,
        results.data(), return_codes.data()), ERROR_VALIDATION__PERCENT_LOW);
}

TEST_F(TestP528Profile, ReturnCodesAreOptional) {
    std::vector<double> d = { 50, 400 };
    std::vector<Result> results(d.size());

    EXPECT_EQ(P528_Profile(d.data(), 2, 10, 1000, 1200, POLARIZATION__HORIZONTAL, 50,
        results.data(), nullptr), SUCCESS);
    EXPECT_EQ(results[0].propagation_mode, PROP_MODE__LOS);
    EXPECT_NE(results[1].propagation_mode, PROP_MODE__LOS);
}

TEST_F(TestP528Profile, WarningsAreSummarized) {
    std::vector<double> d = { 50, 400 };
    std::vector<Result> results(d.size());

    EXPECT_EQ(P528_Profile(d.data(), 2, 3000, 25000, 1200, POLARIZATION__HORIZONTAL, 50,
        results.data(), nullptr), SUCCESS_WITH_WARNINGS);
    EXPECT_TRUE(results[0].warnings & WARNING__HEIGHT_LIMIT_H_2);
}

TEST_F(TestP528Profile, EmptyProfile) {
    EXPECT_EQ(P528_Profile(nullptr, 0, 10, 1000, 1200, POLARIZATION__HORIZONTAL, 50,
        nullptr, nullptr), SUCCESS);
}

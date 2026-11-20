#include "P528GTest.h"
#include "p528.h"
#include "p676.h"
#include "p835.h"

#include <cmath>
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

/*=============================================================================
 |
 |  Description:  Tests for the ray tables behind P528_Profile(). A table
 |                has to reproduce the ray trace it replaces.
 |
 *===========================================================================*/
class TestRayTables : public ::testing::Test {
protected:
    // Allowed deviation of a tabulated ray from the traced one. The traced
    // ray itself is only reproducible to about 1e-8 km, because of rounding
    // in its lowest layers.
    double TOLERANCE__DB = 2.0e-8;
    double TOLERANCE__KM = 2.0e-7;

    RayTraceConfig Config() {
        RayTraceConfig config;
        config.temperature = GlobalTemperature;
        config.dry_pressure = GlobalPressure;
        config.wet_pressure = GlobalWetPressure;
        return config;
    }
};

TEST_F(TestRayTables, SurfaceRayMatchesRayTrace) {
    // 22.235 GHz is the center of a water vapour line, where the absorption
    // is strongest and the kinks of the reference atmosphere show most
    for (double f__ghz : { 0.1, 1.2, 15.5, 22.235, 30.0 }) {
        SpecificAttenuationTable table(f__ghz, Config());
        SurfaceRayTable rays(&table);

        // from 2 m up to the top of the atmosphere, in steps of 0.7 %: this
        // ends rays within every layer and crosses every kink
        int n_tabulated = 0;
        for (double h__km = 0.002; h__km < 99; h__km *= 1.007) {
            SCOPED_TRACE("f__ghz = " + std::to_string(f__ghz) + ", h__km = " + std::to_string(h__km));
            SlantPathAttenuationResult expected;
            SlantPathAttenuation(f__ghz, 0, h__km, PI / 2, &expected, &table);

            double A_gas__db, a__km;
            if (!rays.Lookup(h__km, &A_gas__db, &a__km))
                continue;
            n_tabulated++;
            EXPECT_NEAR(A_gas__db, expected.A_gas__db, TOLERANCE__DB);
            EXPECT_NEAR(a__km, expected.a__km, TOLERANCE__KM);
        }
        // only the first few rays and the very short ones are left to the caller
        EXPECT_GT(n_tabulated, 1100);
    }
}

TEST_F(TestRayTables, LayerRayMatchesRayTrace) {
    struct Geometry { double h_1__km; double h_2__km; double f__ghz; };
    // low and high terminals, and a pair of terminals within one layer
    for (const Geometry& g : { Geometry{ 0.005, 0.03, 0.45 }, Geometry{ 0.01, 1, 1.2 },
                               Geometry{ 1, 12, 22.235 }, Geometry{ 3, 20, 15.5 },
                               Geometry{ 10, 10.02, 30 } }) {
        RayTraceLayers layers;
        GetRayTraceLayers(g.f__ghz, g.h_1__km, g.h_2__km, Config(), &layers);
        LayerRayTable rays(g.h_1__km, &layers);

        // elevation angles from grazing to vertical, in steps of 0.5 %
        int n_tabulated = 0;
        int n_rays = 0;
        for (double theta__rad = 1e-9; theta__rad < PI / 2; theta__rad *= 1.005) {
            SCOPED_TRACE("h_1__km = " + std::to_string(g.h_1__km) + ", theta__rad = " + std::to_string(theta__rad));
            double beta_1__rad = PI / 2 - theta__rad;
            SlantPathAttenuationResult expected;
            RayTraceLoss(g.h_1__km, beta_1__rad, layers, &expected);

            n_rays++;
            double A_gas__db, a__km;
            if (!rays.Lookup(beta_1__rad, &A_gas__db, &a__km))
                continue;
            n_tabulated++;
            EXPECT_NEAR(A_gas__db, expected.A_gas__db, 1e-9 * expected.A_gas__db + 1e-12);
            EXPECT_NEAR(a__km, expected.a__km, 1e-9 * expected.a__km + 1e-12);
        }
        // every range of angles is tabulated after its first few rays
        EXPECT_GT(n_tabulated, 0.8 * n_rays);

        // the end points of the range: a grazing and a vertical ray
        for (double beta_1__rad : { PI / 2, 0.0 }) {
            SlantPathAttenuationResult expected;
            RayTraceLoss(g.h_1__km, beta_1__rad, layers, &expected);
            double A_gas__db, a__km;
            ASSERT_TRUE(rays.Lookup(beta_1__rad, &A_gas__db, &a__km));
            EXPECT_NEAR(A_gas__db, expected.A_gas__db, 1e-9 * expected.A_gas__db + 1e-12);
            EXPECT_NEAR(a__km, expected.a__km, 1e-9 * expected.a__km + 1e-12);
        }
    }
}

TEST_F(TestRayTables, LayerRayLeavesUnsupportedRaysToTheCaller) {
    RayTraceLayers layers;
    GetRayTraceLayers(1.2, 1, 12, Config(), &layers);
    LayerRayTable rays(1, &layers);

    double A_gas__db, a__km;
    // a segment is built once it pays off, not for its first rays
    EXPECT_FALSE(rays.Lookup(1.0, &A_gas__db, &a__km));
    for (int i = 0; i < 30; i++)
        rays.Lookup(1.0 + 1e-3 * i, &A_gas__db, &a__km);
    EXPECT_TRUE(rays.Lookup(1.0, &A_gas__db, &a__km));

    // downward-looking rays do not stay within the layers
    EXPECT_FALSE(rays.Lookup(PI / 2 + 0.01, &A_gas__db, &a__km));
    EXPECT_FALSE(rays.Lookup(-0.01, &A_gas__db, &a__km));
    EXPECT_FALSE(rays.Lookup(std::nan(""), &A_gas__db, &a__km));

    // terminals at the same height: there are no layers to interpolate
    RayTraceLayers no_layers;
    GetRayTraceLayers(1.2, 1, 1, Config(), &no_layers);
    LayerRayTable no_rays(1, &no_layers);
    for (int i = 0; i < 30; i++)
        EXPECT_FALSE(no_rays.Lookup(1.0, &A_gas__db, &a__km));
}

TEST_F(TestRayTables, GrazingRayEqualsSlantPathAttenuation) {
    struct Geometry { double h_1__km; double h_2__km; double f__ghz; double theta_min__rad; };
    // the largest depression angle is the one of the ray that grazes the surface
    for (const Geometry& g : { Geometry{ 0.005, 0.03, 0.45, -0.001 }, Geometry{ 1, 12, 5.1, -0.0145 },
                               Geometry{ 3, 20, 15.5, -0.0248 }, Geometry{ 3, 3, 15.5, -0.0248 } }) {
        SpecificAttenuationTable table(g.f__ghz, Config());
        GrazingRayCache rays(g.f__ghz, g.h_1__km, g.h_2__km, &table);

        // The cache only reorganizes the calculation, so the results are equal.
        // Every angle is asked for twice: first it is traced, then remembered.
        for (int repeat = 0; repeat < 2; repeat++) {
            for (int i = 1; i <= 400; i++) {
                double beta_1__rad = PI / 2 - g.theta_min__rad * i / 400;
                SCOPED_TRACE("h_1__km = " + std::to_string(g.h_1__km) + ", i = " + std::to_string(i));
                SlantPathAttenuationResult expected;
                SlantPathAttenuation(g.f__ghz, g.h_1__km, g.h_2__km, beta_1__rad, &expected, &table);

                double A_gas__db, a__km;
                rays.Lookup(beta_1__rad, &A_gas__db, &a__km);
                EXPECT_EQ(A_gas__db, expected.A_gas__db);
                EXPECT_EQ(a__km, expected.a__km);
            }
        }
    }
}

TEST_F(TestRayTables, SurfaceRayLeavesUnsupportedHeightsToTheCaller) {
    SpecificAttenuationTable table(1.2, Config());
    SurfaceRayTable rays(&table);

    double A_gas__db, a__km;
    // the table is built once it pays off, not for the first rays
    EXPECT_FALSE(rays.Lookup(1.0, &A_gas__db, &a__km));
    for (int i = 0; i < 20; i++)
        rays.Lookup(1.0 + 0.01 * i, &A_gas__db, &a__km);
    EXPECT_TRUE(rays.Lookup(1.0, &A_gas__db, &a__km));

    // heights outside of the reference atmosphere, and rays of a few layers
    EXPECT_FALSE(rays.Lookup(0, &A_gas__db, &a__km));
    EXPECT_FALSE(rays.Lookup(-1, &A_gas__db, &a__km));
    EXPECT_FALSE(rays.Lookup(100, &A_gas__db, &a__km));
    EXPECT_FALSE(rays.Lookup(1e300, &A_gas__db, &a__km));
    EXPECT_FALSE(rays.Lookup(std::nan(""), &A_gas__db, &a__km));
    EXPECT_FALSE(rays.Lookup(1e-4, &A_gas__db, &a__km));
}

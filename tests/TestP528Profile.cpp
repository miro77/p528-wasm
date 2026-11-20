#include "P528GTest.h"
#include "p528.h"
#include "p676.h"
#include "p835.h"
#include "../src/p528/ProfileContext.h"

#include <cmath>
#include <map>
#include <tuple>
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

TEST_F(TestP528Profile, DenseProfileMatchesPointwise) {
    // Tables are only built where enough distances share them. So take many
    // distances from three narrow ranges: upward and downward-looking rays
    // within line-of-sight, and beyond the horizon (which is at 580 km).
    std::vector<double> d;
    for (int i = 0; i < 60; i++) {
        d.push_back(200 + 0.5 * i);
        d.push_back(560 + 0.25 * i);
        d.push_back(700 + 2.0 * i);
    }
    ExpectProfileMatchesPointwise(d, 1000, 12000, 5100, POLARIZATION__HORIZONTAL, 5);
}

TEST_F(TestP528Profile, ParallelMatchesSerial) {
    // enough distances for four threads, with an invalid one in the third range
    std::vector<double> d = MakeDistances(1200, 10000);
    d[6000] = -1;
    std::vector<Result> serial(d.size()), parallel(d.size());
    std::vector<int> serial_rc(d.size()), parallel_rc(d.size());

    int rtn_serial = P528_Profile(d.data(), static_cast<int>(d.size()), 1000, 25000, 5100,
        POLARIZATION__VERTICAL, 5, serial.data(), serial_rc.data());
    int rtn_parallel = P528_ProfileParallel(d.data(), static_cast<int>(d.size()), 1000, 25000, 5100,
        POLARIZATION__VERTICAL, 5, parallel.data(), parallel_rc.data(), 4);

    EXPECT_EQ(rtn_serial, ERROR_VALIDATION__D_KM);
    EXPECT_EQ(rtn_parallel, rtn_serial);
    for (size_t i = 0; i < d.size(); i++) {
        SCOPED_TRACE("d__km = " + std::to_string(d[i]));
        ASSERT_EQ(parallel_rc[i], serial_rc[i]);
        ASSERT_EQ(parallel[i].propagation_mode, serial[i].propagation_mode);
        ASSERT_EQ(parallel[i].warnings, serial[i].warnings);
        // the ranges of the threads build their own tables
        ASSERT_NEAR(parallel[i].A__db, serial[i].A__db, 1e-7);
        ASSERT_NEAR(parallel[i].A_fs__db, serial[i].A_fs__db, 1e-7);
        ASSERT_NEAR(parallel[i].A_a__db, serial[i].A_a__db, 1e-7);
        ASSERT_EQ(parallel[i].d__km, serial[i].d__km);
        ASSERT_EQ(parallel[i].theta_h1__rad, serial[i].theta_h1__rad);
    }

    // return codes stay optional, and small profiles are not split
    EXPECT_EQ(P528_ProfileParallel(d.data(), static_cast<int>(d.size()), 1000, 25000, 5100,
        POLARIZATION__VERTICAL, 5, parallel.data(), nullptr, 0), rtn_serial);
    EXPECT_EQ(P528_ProfileParallel(d.data(), 10, 1000, 25000, 5100,
        POLARIZATION__VERTICAL, 5, parallel.data(), parallel_rc.data(), 4), SUCCESS_WITH_WARNINGS);
    EXPECT_EQ(P528_ProfileParallel(nullptr, 0, 10, 1000, 1200, POLARIZATION__HORIZONTAL, 50,
        nullptr, nullptr, 4), SUCCESS);
}

/*=============================================================================
 |
 |  Description:  P528_Profile() against the data tables of Recommendation
 |                ITU-R P.528-5. Each column of a table is a profile: the
 |                loss versus distance for one pair of terminal heights.
 |                Unlike the tests of P528(), which sample the tables, these
 |                compare every value.
 |
 *===========================================================================*/
class TestP528ProfileDataTables : public ::testing::TestWithParam<std::string> {
protected:
    // the tables are rounded to 0.1 dB
    double TOLERANCE = 1.0e-1;
};

TEST_P(TestP528ProfileDataTables, EveryColumnIsAProfile) {
    std::vector<InputsAndResult> table = ReadP528InputsAndResultFromDataTable(GetParam(), 1);
    ASSERT_NE(table.size(), 0u);

    // collect the columns: rows of equal heights, in the order of the table
    typedef std::tuple<double, double, double, int, double> Column;
    std::map<Column, std::vector<size_t>> columns;
    for (size_t i = 0; i < table.size(); i++)
        columns[Column(table[i].h_1__meter, table[i].h_2__meter, table[i].f__mhz,
            table[i].T_pol, table[i].p)].push_back(i);

    for (const auto& column : columns) {
        const std::vector<size_t>& rows = column.second;
        const InputsAndResult& inputs = table[rows[0]];
        std::vector<double> d(rows.size());
        for (size_t k = 0; k < rows.size(); k++)
            d[k] = table[rows[k]].d__km;

        std::vector<Result> results(d.size());
        std::vector<int> return_codes(d.size());
        EXPECT_EQ(P528_Profile(d.data(), static_cast<int>(d.size()), inputs.h_1__meter,
            inputs.h_2__meter, inputs.f__mhz, inputs.T_pol, inputs.p, results.data(),
            return_codes.data()), SUCCESS);

        for (size_t k = 0; k < rows.size(); k++) {
            SCOPED_TRACE(testing::Message() << "h_1__meter=" << inputs.h_1__meter << ", h_2__meter="
                << inputs.h_2__meter << ", f__mhz=" << inputs.f__mhz << ", p=" << inputs.p
                << ", d__km=" << d[k]);
            EXPECT_EQ(return_codes[k], SUCCESS);
            EXPECT_NEAR(results[k].A__db, table[rows[k]].expectedResult.A__db, TOLERANCE);
            if (!std::isnan(table[rows[k]].expectedResult.A_fs__db))
                EXPECT_NEAR(results[k].A_fs__db, table[rows[k]].expectedResult.A_fs__db, TOLERANCE);
        }
    }
}

INSTANTIATE_TEST_SUITE_P(TestP528, TestP528ProfileDataTables,
    ::testing::ValuesIn(GetDataTablesFileList()));

/*=============================================================================
 |
 |  Description:  Tests for the search of the reflection angle that
 |                P528_Profile() uses. It has to end at exactly the angle
 |                FindPsiAtDistance() finds: the search stops within 1 m of
 |                the distance, so another angle would be another result.
 |
 *===========================================================================*/
class TestPsiSearch : public ::testing::Test {
protected:
    void SetTerminals(double h_1__meter, double h_2__meter, double f__mhz) {
        terminal_1.h_r__km = h_1__meter / 1000;
        terminal_2.h_r__km = h_2__meter / 1000;
        TerminalGeometry(f__mhz, &terminal_1);
        TerminalGeometry(f__mhz, &terminal_2);
    }

    Terminal terminal_1;
    Terminal terminal_2;
};

TEST_F(TestPsiSearch, DistanceEqualsRayOptics) {
    for (double h_2__meter : { 30.0, 1000.0, 20000.0 }) {
        SetTerminals(5, h_2__meter, 1200);
        // angles from grazing to vertical, dense near grazing
        for (int i = 0; i < 20000; i++) {
            double psi = (PI / 2) * pow((i + 0.5) / 20000, 3);
            LineOfSightParams params;
            RayOptics(&terminal_1, &terminal_2, psi, &params);
            ASSERT_EQ(RayOpticsDistance(&terminal_1, &terminal_2, psi), params.d__km)
                << "psi = " << psi;
        }
    }
}

TEST_F(TestPsiSearch, FindsTheAngleOfFindPsiAtDistance) {
    SetTerminals(1000, 12000, 5100);
    double d_ML__km = terminal_1.d_r__km + terminal_2.d_r__km;

    PsiSearchTree tree;
    // distances in an arbitrary order, twice: the second pass only reads the tree
    for (int repeat = 0; repeat < 2; repeat++) {
        for (int i = 0; i < 3000; i++) {
            double d__km = d_ML__km * ((i * 7919) % 3000) / 3000.0;
            ASSERT_EQ(FindPsiAtDistance(d__km, &terminal_1, &terminal_2, &tree),
                FindPsiAtDistance(d__km, &terminal_1, &terminal_2)) << "d__km = " << d__km;
        }
    }
    // beyond the horizon the search does not converge, it ends at its smallest step
    EXPECT_EQ(FindPsiAtDistance(d_ML__km + 10, &terminal_1, &terminal_2, &tree),
        FindPsiAtDistance(d_ML__km + 10, &terminal_1, &terminal_2));
    EXPECT_LE(tree.nodes.size(), PsiSearchTree::MAX_NODES);
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

TEST_F(TestRayTables, TerminalGeometryOnlyInterpolatesTheAbsorption) {
    for (double f__mhz : { 300.0, 5100.0, 22235.0 }) {
        SpecificAttenuationTable table(f__mhz / 1000, Config());
        for (double h_r__km : { 0.0015, 0.03, 1.0, 12.0, 20.0, 80.0 }) {
            SCOPED_TRACE("f__mhz = " + std::to_string(f__mhz) + ", h_r__km = " + std::to_string(h_r__km));
            Terminal exact, tabulated;
            exact.h_r__km = h_r__km;
            tabulated.h_r__km = h_r__km;
            TerminalGeometry(f__mhz, &exact);
            TerminalGeometry(f__mhz, &tabulated, &table);

            // the geometry decides propagation modes and search steps: no rounding allowed
            EXPECT_EQ(tabulated.d_r__km, exact.d_r__km);
            EXPECT_EQ(tabulated.a__km, exact.a__km);
            EXPECT_EQ(tabulated.theta__rad, exact.theta__rad);
            EXPECT_EQ(tabulated.phi__rad, exact.phi__rad);
            EXPECT_EQ(tabulated.h_e__km, exact.h_e__km);
            EXPECT_EQ(tabulated.delta_h__km, exact.delta_h__km);

            EXPECT_NEAR(tabulated.A_a__db, exact.A_a__db, 1e-9 * exact.A_a__db + 1e-12);
        }
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

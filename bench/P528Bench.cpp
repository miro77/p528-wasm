/** @file P528Bench.cpp
 * Benchmark for evaluating P.528 along a profile (many distances that share
 * the same terminal heights, frequency, polarization and time percentage).
 *
 * For every scenario the benchmark
 *   1. times the point-by-point reference (one P528() call per distance),
 *   2. times the profile API (one P528_Profile() call for all distances),
 *   3. compares both result sets, so a speed-up never silently changes results.
 *
 * The point-by-point results can additionally be saved to / checked against a
 * reference file. That guards optimizations which also change the
 * point-by-point code path (where comparing the two APIs proves nothing).
 *
 * Usage: P528Bench [--points N] [--repeat R] [--csv FILE] [--stages]
 *                  [--profile-only]
 *                  [--save-reference FILE] [--check-reference FILE]
 */
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

#include "p528.h"
#include "p676.h"
#include "p835.h"

namespace ITS {
namespace ITU {
namespace PSeries {
namespace p528 {
// internal function of LineOfSight.cpp, only needed for the stage breakdown
double FindPsiAtDistance(double d__km, Terminal* terminal_1,
                         Terminal* terminal_2);
}  // namespace p528
}  // namespace PSeries
}  // namespace ITU
}  // namespace ITS

using namespace ITS::ITU::PSeries::p528;

namespace {

using Clock = std::chrono::steady_clock;

double SecondsSince(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}

struct Scenario {
  const char* name;
  double h_1__meter;
  double h_2__meter;
  double f__mhz;
  int T_pol;
  double p;
  double d_max__km;
};

// The scenarios cover the three propagation regions (line-of-sight,
// diffraction, troposcatter) for low, medium and high terminals.
const std::vector<Scenario> kScenarios = {
    {"ground-air 10m/1km 1.2GHz", 10, 1000, 1200, POLARIZATION__HORIZONTAL, 50,
     600},
    {"ground-air 15m/10km 300MHz", 15, 10000, 300, POLARIZATION__VERTICAL, 10,
     1000},
    {"air-air 1km/12km 5.1GHz", 1000, 12000, 5100, POLARIZATION__HORIZONTAL, 5,
     1200},
    {"low-low 5m/30m 450MHz", 5, 30, 450, POLARIZATION__VERTICAL, 95, 300},
    {"high 3km/20km 15.5GHz", 3000, 20000, 15500, POLARIZATION__VERTICAL, 50,
     1800},
};

std::vector<double> MakeProfile(double d_max__km, int n) {
  std::vector<double> d(n);
  for (int i = 0; i < n; i++)
    d[i] = (n == 1) ? d_max__km : d_max__km * i / (n - 1);
  return d;
}

struct Comparison {
  double max_abs_diff__db = 0;   // largest |A__db| deviation
  double max_abs_diff_d__km = 0; // distance at which it occurs
  int rc_mismatches = 0;         // differing return codes
  int mode_mismatches = 0;       // differing propagation modes
  int warning_mismatches = 0;    // differing warning flags
};

Comparison Compare(const std::vector<double>& d,
                   const std::vector<Result>& ref,
                   const std::vector<int>& ref_rc,
                   const std::vector<Result>& test,
                   const std::vector<int>& test_rc) {
  Comparison c;
  for (size_t i = 0; i < d.size(); i++) {
    if (ref_rc[i] != test_rc[i]) c.rc_mismatches++;
    if (ref[i].propagation_mode != test[i].propagation_mode)
      c.mode_mismatches++;
    if (ref[i].warnings != test[i].warnings) c.warning_mismatches++;
    const double diffs[] = {
        std::fabs(ref[i].A__db - test[i].A__db),
        std::fabs(ref[i].A_fs__db - test[i].A_fs__db),
        std::fabs(ref[i].A_a__db - test[i].A_a__db),
    };
    for (double diff : diffs) {
      if (diff > c.max_abs_diff__db || std::isnan(diff)) {
        c.max_abs_diff__db = diff;
        c.max_abs_diff_d__km = d[i];
      }
    }
  }
  return c;
}

/** One reference line per profile point; %.17g round-trips a double. */
void SaveReference(std::FILE* file, int scenario, const std::vector<double>& d,
                   const std::vector<Result>& res, const std::vector<int>& rc) {
  for (size_t i = 0; i < d.size(); i++)
    std::fprintf(file, "%d %.17g %d %d %d %.17g %.17g %.17g %.17g\n", scenario,
                 d[i], rc[i], res[i].propagation_mode, res[i].warnings,
                 res[i].A__db, res[i].A_fs__db, res[i].A_a__db,
                 res[i].theta_h1__rad);
}

/** Reads the reference results of one scenario; false on a format mismatch. */
bool LoadReference(std::FILE* file, int scenario, const std::vector<double>& d,
                   std::vector<Result>* res, std::vector<int>* rc) {
  for (size_t i = 0; i < d.size(); i++) {
    int s_idx;
    double d_ref;
    Result& r = (*res)[i];
    if (std::fscanf(file, "%d %lf %d %d %d %lf %lf %lf %lf", &s_idx, &d_ref,
                    &(*rc)[i], &r.propagation_mode, &r.warnings, &r.A__db,
                    &r.A_fs__db, &r.A_a__db, &r.theta_h1__rad) != 9 ||
        s_idx != scenario || std::fabs(d_ref - d[i]) > 1e-9)
      return false;
    r.d__km = d_ref;
  }
  return true;
}

/** Times the building blocks of a single P528() call for one scenario. */
void PrintStageBreakdown(const Scenario& s, int repeat) {
  auto time_it = [repeat](const char* label, const std::function<void()>& fn) {
    auto start = Clock::now();
    for (int i = 0; i < repeat; i++) fn();
    std::printf("    %-46s %10.1f us\n", label,
                1e6 * SecondsSince(start) / repeat);
  };

  Terminal t1, t2;
  t1.h_r__km = s.h_1__meter / 1000;
  t2.h_r__km = s.h_2__meter / 1000;
  std::printf("  stage breakdown (per call):\n");
  time_it("TerminalGeometry (low + high terminal)", [&] {
    TerminalGeometry(s.f__mhz, &t1);
    TerminalGeometry(s.f__mhz, &t2);
  });
  const double d_ML__km = t1.d_r__km + t2.d_r__km;

  SlantPathAttenuationResult slant;
  time_it("SlantPathAttenuation (h_1 -> h_2, one ray)", [&] {
    SlantPathAttenuation(s.f__mhz / 1000, t1.h_r__km, t2.h_r__km, PI / 2 - 0.01,
                         &slant);
  });
  time_it("SmoothEarthDiffraction (x2, diffraction line)", [&] {
    SmoothEarthDiffraction(t1.d_r__km, t2.d_r__km, s.f__mhz, d_ML__km + 10,
                           s.T_pol);
    SmoothEarthDiffraction(t1.d_r__km, t2.d_r__km, s.f__mhz, d_ML__km + 30,
                           s.T_pol);
  });
  TroposcatterParams tropo;
  time_it("Troposcatter (one distance)", [&] {
    Troposcatter(&t1, &t2, d_ML__km + 50, s.f__mhz, &tropo);
  });
  Result r;
  time_it("P528 full call, line-of-sight point", [&] {
    P528(0.5 * d_ML__km, s.h_1__meter, s.h_2__meter, s.f__mhz, s.T_pol, s.p,
         &r);
  });
  time_it("P528 full call, transhorizon point", [&] {
    P528(d_ML__km + 50, s.h_1__meter, s.h_2__meter, s.f__mhz, s.T_pol, s.p, &r);
  });

  // Building blocks that remain per point in P528_Profile()
  std::printf("  stage breakdown of the profile path (per point):\n");
  const int fast_repeat = 200 * repeat;
  auto time_fast = [fast_repeat](const char* label,
                                 const std::function<void(int)>& fn) {
    auto start = Clock::now();
    for (int i = 0; i < fast_repeat; i++) fn(i);
    std::printf("    %-46s %10.2f us\n", label,
                1e6 * SecondsSince(start) / fast_repeat);
  };
  RayTraceConfig config;
  config.temperature = GlobalTemperature;
  config.dry_pressure = GlobalPressure;
  config.wet_pressure = GlobalWetPressure;
  SpecificAttenuationTable table(s.f__mhz / 1000, config);
  RayTraceLayers layers;
  GetRayTraceLayers(s.f__mhz / 1000, t1.h_r__km, t2.h_r__km, config, &layers);
  volatile double sink = 0;
  time_fast("LOS: FindPsiAtDistance", [&](int i) {
    sink = FindPsiAtDistance(d_ML__km * (0.2 + 0.6 * i / fast_repeat), &t1, &t2);
  });
  time_fast("LOS: ray through shared layers", [&](int i) {
    RayTrace(t1.h_r__km, PI / 2 - 0.001 - 0.2 * i / fast_repeat, layers,
             &slant);
  });
  time_fast("LOS: grazing ray (table)", [&](int i) {
    SlantPathAttenuation(s.f__mhz / 1000, t1.h_r__km, t2.h_r__km,
                         PI / 2 + 1e-4 * (1 + i % 7), &slant, &table);
  });
  time_fast("transhorizon: Troposcatter", [&](int i) {
    Troposcatter(&t1, &t2, d_ML__km + 50 + 0.01 * i, s.f__mhz, &tropo);
  });
  time_fast("transhorizon: ray 0 -> h_v (table)", [&](int i) {
    SlantPathAttenuation(s.f__mhz / 1000, 0, tropo.h_v__km * (1 + 1e-4 * i),
                         PI / 2, &slant, &table);
  });
  time_fast("variability: LongTermVariability (x2)", [&](int i) {
    double Y, A_Y;
    LongTermVariability(t1.d_r__km, t2.d_r__km, d_ML__km + 0.01 * i, s.f__mhz,
                        s.p, 1, -20, &Y, &A_Y);
    LongTermVariability(t1.d_r__km, t2.d_r__km, d_ML__km + 0.01 * i, s.f__mhz,
                        50, 1, -20, &Y, &A_Y);
  });
  time_fast("variability: NakagamiRice + Combine", [&](int i) {
    sink = CombineDistributions(1, 3, 0, NakagamiRice(5 + 1e-3 * i, s.p), s.p);
  });
  (void)sink;
}

}  // namespace

int main(int argc, char** argv) {
  int n_points = 2000;
  int repeat = 1;
  bool stages = false;
  bool profile_only = false;
  std::string csv_path, save_path, check_path;

  for (int i = 1; i < argc; i++) {
    if (!std::strcmp(argv[i], "--points") && i + 1 < argc)
      n_points = std::atoi(argv[++i]);
    else if (!std::strcmp(argv[i], "--repeat") && i + 1 < argc)
      repeat = std::atoi(argv[++i]);
    else if (!std::strcmp(argv[i], "--csv") && i + 1 < argc)
      csv_path = argv[++i];
    else if (!std::strcmp(argv[i], "--save-reference") && i + 1 < argc)
      save_path = argv[++i];
    else if (!std::strcmp(argv[i], "--check-reference") && i + 1 < argc)
      check_path = argv[++i];
    else if (!std::strcmp(argv[i], "--profile-only"))
      profile_only = true;
    else if (!std::strcmp(argv[i], "--stages"))
      stages = true;
    else {
      std::printf(
          "Usage: %s [--points N] [--repeat R] [--csv FILE] [--stages]\n"
          "          [--profile-only]\n"
          "          [--save-reference FILE] [--check-reference FILE]\n"
          "  --points N  number of profile points per scenario (default 2000)\n"
          "  --repeat R  repetitions, the fastest run is reported (default 1)\n"
          "  --csv FILE  write the per-scenario timings to FILE\n"
          "  --stages    also time the building blocks of a single call\n"
          "  --profile-only  skip the slow point-by-point run, for timing\n"
          "                  large profiles (no comparison of the two APIs)\n"
          "  --save-reference FILE   store the point-by-point results\n"
          "  --check-reference FILE  compare both APIs against stored results\n"
          "                          (needs the same --points as when saved)\n",
          argv[0]);
      return !std::strcmp(argv[i], "--help") ? 0 : 1;
    }
  }
  if (n_points < 1 || repeat < 1) {
    std::printf("--points and --repeat must be >= 1\n");
    return 1;
  }
  if (profile_only && (!save_path.empty() || !check_path.empty())) {
    std::printf("--profile-only cannot be combined with a reference file\n");
    return 1;
  }

  std::printf("P.528 profile benchmark: %d points per scenario, best of %d\n\n",
              n_points, repeat);
  std::printf("%-28s %12s %12s %9s %12s %s\n", "scenario", "point [ms]",
              "profile [ms]", "speed-up", "max diff[dB]", "mismatches");

  std::ofstream csv;
  if (!csv_path.empty()) {
    csv.open(csv_path);
    csv << "scenario,points,pointwise_ms,profile_ms,speedup,max_abs_diff_db,"
           "max_diff_d_km,rc_mismatches,mode_mismatches,warning_mismatches\n";
  }

  std::FILE* save_file = nullptr;
  std::FILE* check_file = nullptr;
  if (!save_path.empty() && !(save_file = std::fopen(save_path.c_str(), "w"))) {
    std::printf("cannot write %s\n", save_path.c_str());
    return 1;
  }
  if (!check_path.empty() &&
      !(check_file = std::fopen(check_path.c_str(), "r"))) {
    std::printf("cannot read %s\n", check_path.c_str());
    return 1;
  }

  double total_ref = 0, total_profile = 0, worst_diff = 0;
  double worst_stored_diff = 0;
  int total_mismatches = 0;

  for (size_t i_scenario = 0; i_scenario < kScenarios.size(); i_scenario++) {
    const Scenario& s = kScenarios[i_scenario];
    const std::vector<double> d = MakeProfile(s.d_max__km, n_points);
    std::vector<Result> ref(n_points), prof(n_points);
    std::vector<int> ref_rc(n_points), prof_rc(n_points);

    double t_ref = 1e300, t_prof = 1e300;
    for (int r = 0; r < repeat; r++) {
      auto start = Clock::now();
      if (!profile_only) {
        for (int i = 0; i < n_points; i++)
          ref_rc[i] = P528(d[i], s.h_1__meter, s.h_2__meter, s.f__mhz, s.T_pol,
                           s.p, &ref[i]);
        t_ref = std::fmin(t_ref, SecondsSince(start));
      }

      start = Clock::now();
      P528_Profile(d.data(), n_points, s.h_1__meter, s.h_2__meter, s.f__mhz,
                   s.T_pol, s.p, prof.data(), prof_rc.data());
      t_prof = std::fmin(t_prof, SecondsSince(start));
    }

    if (profile_only) {
      std::printf("%-28s %12s %12.1f %9s %12s %s   (%.1f us/point)\n", s.name,
                  "-", 1e3 * t_prof, "-", "-", "-", 1e6 * t_prof / n_points);
      total_profile += t_prof;
      continue;
    }

    const Comparison c = Compare(d, ref, ref_rc, prof, prof_rc);
    const int mismatches =
        c.rc_mismatches + c.mode_mismatches + c.warning_mismatches;
    std::printf("%-28s %12.1f %12.1f %8.1fx %12.3g %d\n", s.name, 1e3 * t_ref,
                1e3 * t_prof, t_ref / t_prof, c.max_abs_diff__db, mismatches);
    if (csv.is_open())
      csv << s.name << ',' << n_points << ',' << 1e3 * t_ref << ','
          << 1e3 * t_prof << ',' << t_ref / t_prof << ',' << c.max_abs_diff__db
          << ',' << c.max_abs_diff_d__km << ',' << c.rc_mismatches << ','
          << c.mode_mismatches << ',' << c.warning_mismatches << '\n';

    if (save_file) SaveReference(save_file, (int)i_scenario, d, ref, ref_rc);
    if (check_file) {
      std::vector<Result> stored(n_points);
      std::vector<int> stored_rc(n_points);
      if (!LoadReference(check_file, (int)i_scenario, d, &stored, &stored_rc)) {
        std::printf("reference file %s does not match this run (different "
                    "--points or scenarios?)\n", check_path.c_str());
        return 1;
      }
      const Comparison c_point = Compare(d, stored, stored_rc, ref, ref_rc);
      const Comparison c_prof = Compare(d, stored, stored_rc, prof, prof_rc);
      std::printf("  vs stored reference: point-by-point %.3g dB (at %.3f km),"
                  " profile %.3g dB (at %.3f km)\n",
                  c_point.max_abs_diff__db, c_point.max_abs_diff_d__km,
                  c_prof.max_abs_diff__db, c_prof.max_abs_diff_d__km);
      for (const Comparison* cmp : {&c_point, &c_prof}) {
        total_mismatches += cmp->rc_mismatches + cmp->mode_mismatches +
                            cmp->warning_mismatches;
        if (cmp->max_abs_diff__db > worst_stored_diff ||
            std::isnan(cmp->max_abs_diff__db))
          worst_stored_diff = cmp->max_abs_diff__db;
      }
    }

    if (stages) PrintStageBreakdown(s, 20);

    total_ref += t_ref;
    total_profile += t_prof;
    total_mismatches += mismatches;
    if (c.max_abs_diff__db > worst_diff || std::isnan(c.max_abs_diff__db))
      worst_diff = c.max_abs_diff__db;
  }

  if (profile_only) {
    std::printf("%-28s %12s %12.1f\n", "TOTAL", "-", 1e3 * total_profile);
    return 0;
  }
  std::printf("%-28s %12.1f %12.1f %8.1fx %12.3g %d\n", "TOTAL",
              1e3 * total_ref, 1e3 * total_profile, total_ref / total_profile,
              worst_diff, total_mismatches);

  if (save_file) std::fclose(save_file);
  if (check_file) {
    std::fclose(check_file);
    std::printf("largest deviation from the stored reference: %.3g dB\n",
                worst_stored_diff);
  }

  // A non-zero exit code flags a result regression to scripts and CI.
  const double kTolerance__db = 1e-6;
  if (total_mismatches != 0 || !(worst_diff <= kTolerance__db) ||
      !(worst_stored_diff <= kTolerance__db)) {
    std::printf("\nFAILED: results deviate from the reference "
                "(tolerance %g dB)\n", kTolerance__db);
    return 2;
  }
  return 0;
}

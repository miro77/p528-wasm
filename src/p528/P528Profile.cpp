/** @file P528Profile.cpp
 * Evaluates the P.528 model for many distances along one profile.
 */
#include "p528.h"

#include <vector>

#include "ProfileContext.h"
#include "p835.h"

// WebAssembly has no threads unless the page is served with the headers that
// enable SharedArrayBuffer. The module is built without them.
#if defined(__EMSCRIPTEN__) && !defined(__EMSCRIPTEN_PTHREADS__)
#define P528_SINGLE_THREADED
#else
#include <thread>
#endif

namespace ITS {
namespace ITU {
namespace PSeries {
namespace p528 {

namespace {

// A worker thread builds its own tables, which costs about as much as a few
// hundred distances. Below this number of distances per thread it is not
// worth starting one.
const int MIN_POINTS_PER_THREAD = 2048;

/** Adds the return code of a distance (or range) to the summary of a profile:
 *  errors take precedence over warnings, the first error is kept. */
void Summarize(int rc, int* summary) {
  if (rc != SUCCESS &&
      (*summary == SUCCESS ||
       (*summary == SUCCESS_WITH_WARNINGS && rc != SUCCESS_WITH_WARNINGS)))
    *summary = rc;
}

/**
 * Evaluates a range of distances one after the other. The state of the
 * profile (terminal geometry, tables) belongs to the call, so ranges can be
 * evaluated at the same time.
 *
 * The parameters are those of P528_Profile().
 ****************************************************************************/
int EvaluateRange(const double* d__km, int n_points, double h_1__meter,
                  double h_2__meter, double f__mhz, int T_pol, double p,
                  Result* results, int* return_codes) {
  ProfileContext context;
  bool context_ready = false;
  int summary = SUCCESS;

  for (int i = 0; i < n_points; i++) {
    Result* result = &results[i];
    result->A_fs__db = 0;
    result->A_a__db = 0;
    result->A__db = 0;
    result->d__km = 0;
    result->theta_h1__rad = 0;
    result->propagation_mode = PROP_MODE__NOT_SET;
    result->warnings = WARNING__NO_WARNINGS;

    // Validation depends on the distance too, so it runs for every point.
    int rc = ValidateInputs(d__km[i], h_1__meter, h_2__meter, f__mhz, T_pol, p,
                            &result->warnings);
    if (rc == ERROR_HEIGHT_AND_DISTANCE) {
      rc = SUCCESS;  // same special case as in P528_Ex()
    } else if (rc == SUCCESS) {
      // Reaching this point means that the distance-independent inputs are
      // valid, so the shared state can be built (once).
      if (!context_ready) {
        // All rays of the profile share one attenuation table for the
        // frequency, including those of the terminal geometry.
        RayTraceConfig config;
        config.temperature = GlobalTemperature;
        config.dry_pressure = GlobalPressure;
        config.wet_pressure = GlobalWetPressure;
        context.gamma_table.reset(
            new SpecificAttenuationTable(f__mhz / 1000, config));
        context.surface_rays.reset(
            new SurfaceRayTable(context.gamma_table.get()));

        InitializeProfileContext(h_1__meter, h_2__meter, f__mhz, T_pol, p,
                                 &context);
        context_ready = true;
      }

      TroposcatterParams tropo;
      LineOfSightParams los_params;
      EvaluateProfilePoint(&context, d__km[i], result, &tropo, &los_params);

      if (result->warnings != WARNING__NO_WARNINGS) rc = SUCCESS_WITH_WARNINGS;
    }

    if (return_codes) return_codes[i] = rc;

    Summarize(rc, &summary);
  }

  return summary;
}

}  // namespace

/**
 * Calculates the basic transmission loss for many path distances that share
 * the terminal heights, frequency, polarization and time percentage, e.g. the
 * points of a profile.
 *
 * Every point gets the result P528() would return for it, within 1e-6 dB. The
 * function is much faster than calling P528() per point:
 *   - Everything that does not depend on the distance is computed only once.
 *   - The specific gaseous attenuation is tabulated versus height instead of
 *     being summed line by line for every layer of every ray.
 *   - The rays themselves are tabulated where their loss is a smooth function
 *     of one variable, and remembered where it is not.
 *
 * The function runs on the calling thread. See P528_ProfileParallel() for
 * large profiles.
 *
 * @param[in]  d__km         %Path distances, in km (n_points values, any order)
 * @param[in]  n_points      Number of distances
 * @param[in]  h_1__meter    Height of the low terminal, in meters
 * @param[in]  h_2__meter    Height of the high terminal, in meters
 * @param[in]  f__mhz        Frequency, in MHz
 * @param[in]  T_pol         Polarization type (POLARIZATION__HORIZONTAL or
 *                           POLARIZATION__VERTICAL)
 * @param[in]  p             Time percentage
 * @param[out] results       Result structures (n_points values)
 * @param[out] return_codes  RETURNCODES value per distance (n_points values).
 *                           May be null if only the summary is of interest.
 * @return                   SUCCESS if every point succeeded. Otherwise the
 *                           first error code, or SUCCESS_WITH_WARNINGS if
 *                           there were only warnings.
 ****************************************************************************/
int P528_Profile(const double* d__km, int n_points, double h_1__meter,
                 double h_2__meter, double f__mhz, int T_pol, double p,
                 Result* results, int* return_codes) {
  return EvaluateRange(d__km, n_points, h_1__meter, h_2__meter, f__mhz, T_pol,
                       p, results, return_codes);
}

/**
 * P528_Profile() on several threads, for large profiles.
 *
 * The distances are split into contiguous ranges, one per thread. Every range
 * is a profile of its own, with its own state, so the threads share nothing
 * but the input and output arrays. Each thread therefore pays the fixed cost
 * of a profile, and a range has to be large enough to make up for that:
 * profiles of less than 4096 distances are evaluated on the calling thread.
 *
 * The results are those of P528_Profile(), up to the accuracy of its tables
 * (about 1e-8 dB): a table is built from the distances that use it, which
 * differ from range to range.
 *
 * In a WebAssembly build the function is P528_Profile().
 *
 * @param[in]  d__km         %Path distances, in km (n_points values, any order)
 * @param[in]  n_points      Number of distances
 * @param[in]  h_1__meter    Height of the low terminal, in meters
 * @param[in]  h_2__meter    Height of the high terminal, in meters
 * @param[in]  f__mhz        Frequency, in MHz
 * @param[in]  T_pol         Polarization type (POLARIZATION__HORIZONTAL or
 *                           POLARIZATION__VERTICAL)
 * @param[in]  p             Time percentage
 * @param[out] results       Result structures (n_points values)
 * @param[out] return_codes  RETURNCODES value per distance (n_points values).
 *                           May be null if only the summary is of interest.
 * @param[in]  max_threads   Largest number of threads to use, including the
 *                           calling one. 0 or less: one per processor.
 * @return                   SUCCESS if every point succeeded. Otherwise the
 *                           first error code, or SUCCESS_WITH_WARNINGS if
 *                           there were only warnings.
 ****************************************************************************/
int P528_ProfileParallel(const double* d__km, int n_points, double h_1__meter,
                         double h_2__meter, double f__mhz, int T_pol, double p,
                         Result* results, int* return_codes, int max_threads) {
#ifdef P528_SINGLE_THREADED
  (void)max_threads;
  return EvaluateRange(d__km, n_points, h_1__meter, h_2__meter, f__mhz, T_pol,
                       p, results, return_codes);
#else
  if (max_threads <= 0) {
    // 0 means that the number of processors is not known
    max_threads = (int)std::thread::hardware_concurrency();
    if (max_threads <= 0) max_threads = 1;
  }
  const int n_ranges = MIN(max_threads, n_points / MIN_POINTS_PER_THREAD);
  if (n_ranges <= 1)
    return EvaluateRange(d__km, n_points, h_1__meter, h_2__meter, f__mhz,
                         T_pol, p, results, return_codes);

  // range r covers the distances first[r] ... first[r + 1] - 1
  std::vector<int> first(n_ranges + 1);
  for (int r = 0; r <= n_ranges; r++)
    first[r] = (int)((long long)n_points * r / n_ranges);

  std::vector<int> summaries(n_ranges, SUCCESS);
  auto evaluate = [&](int r) {
    summaries[r] = EvaluateRange(
        d__km + first[r], first[r + 1] - first[r], h_1__meter, h_2__meter,
        f__mhz, T_pol, p, results + first[r],
        return_codes ? return_codes + first[r] : nullptr);
  };

  // The calling thread takes the first range. If the system refuses to start
  // a thread, its range is evaluated on the calling thread as well.
  std::vector<std::thread> threads;
  std::vector<int> not_started;
  threads.reserve(n_ranges - 1);
  for (int r = 1; r < n_ranges; r++) {
    try {
      threads.emplace_back(evaluate, r);
    } catch (...) {
      not_started.push_back(r);
    }
  }
  evaluate(0);
  for (int r : not_started) evaluate(r);
  for (std::thread& thread : threads) thread.join();

  // the first error of the first range with an error is the first of all
  int summary = SUCCESS;
  for (int r = 0; r < n_ranges; r++) Summarize(summaries[r], &summary);
  return summary;
#endif
}

}  // namespace p528
}  // namespace PSeries
}  // namespace ITU
}  // namespace ITS

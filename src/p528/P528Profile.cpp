/** @file P528Profile.cpp
 * Evaluates the P.528 model for many distances along one profile.
 */
#include "p528.h"

#include "ProfileContext.h"
#include "p835.h"

namespace ITS {
namespace ITU {
namespace PSeries {
namespace p528 {

/**
 * Calculates the basic transmission loss for many path distances that share
 * the terminal heights, frequency, polarization and time percentage, e.g. the
 * points of a profile.
 *
 * Every point gets the result P528() would return for it. The function is
 * much faster than calling P528() per point, because everything that does not
 * depend on the distance is computed only once, and because the specific
 * gaseous attenuation is tabulated versus height instead of being summed line
 * by line for every ray. The interpolation changes the atmospheric absorption
 * by less than 1e-6 dB.
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
        InitializeProfileContext(h_1__meter, h_2__meter, f__mhz, T_pol, p,
                                 &context);
        context_ready = true;

        // All quantities above are exact. The rays traced per distance share
        // one attenuation table for the frequency.
        RayTraceConfig config;
        config.temperature = GlobalTemperature;
        config.dry_pressure = GlobalPressure;
        config.wet_pressure = GlobalWetPressure;
        context.gamma_table.reset(
            new SpecificAttenuationTable(f__mhz / 1000, config));
      }

      TroposcatterParams tropo;
      LineOfSightParams los_params;
      EvaluateProfilePoint(&context, d__km[i], result, &tropo, &los_params);

      if (result->warnings != WARNING__NO_WARNINGS) rc = SUCCESS_WITH_WARNINGS;
    }

    if (return_codes) return_codes[i] = rc;

    // errors take precedence over warnings in the summary
    if (rc != SUCCESS &&
        (summary == SUCCESS || (summary == SUCCESS_WITH_WARNINGS &&
                                rc != SUCCESS_WITH_WARNINGS)))
      summary = rc;
  }

  return summary;
}

}  // namespace p528
}  // namespace PSeries
}  // namespace ITU
}  // namespace ITS

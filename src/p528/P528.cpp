/** @file P528.cpp
 *  Implements the model from ITU-R P.528-5. ANnex 2, Section 3.
 */

#include "p528.h"

#include <math.h>

#include "p676.h"
#include "ProfileContext.h"

namespace ITS {
namespace ITU {
namespace PSeries {
namespace p528 {

/**
 * This is the main entry point to this software. It calculates the basic
 * transmission loss for a user-specified path using Recommendation ITU-R
 * P.528-5
 *
 * References:
 *     - Recommendation ITU‑R P.528‑5 "Propagation curves for aeronautical
 * mobile and radionavigation services using the VHF, UHF and SHF bands", Annex
 * 2, Section 3.
 *
 * @param[in]  d__km       %Path distance, in km
 * @param[in]  h_1__meter  Height of the low terminal, in meters
 * @param[in]  h_2__meter  Height of the high terminal, in meters
 * @param[in]  f__mhz      Frequency, in MHz
 * @param[in]  T_pol       Polarization type (POLARIZATION__HORIZONTAL or
 *                         POLARIZATION__VERTICAL)
 * @param[in]  p           Time percentage
 * @param[out] result      Result structure containing various
 * @return                 RETURNCODES value explaining success of code run
 ****************************************************************************/
int P528(double d__km, double h_1__meter, double h_2__meter, double f__mhz,
         int T_pol, double p, Result* result) {
  Terminal terminal_1;
  Terminal terminal_2;
  TroposcatterParams tropo;
  Path path;
  LineOfSightParams los_params;

  return P528_Ex(d__km, h_1__meter, h_2__meter, f__mhz, T_pol, p, result,
                 &terminal_1, &terminal_2, &tropo, &path, &los_params);
}

/**
 * This is the main entry point to this software. It calculates the basic
 * transmission loss for a user-specified path using Recommendation ITU-R
 * P.528-5
 *
 * References:
 *     - Recommendation ITU‑R P.528‑5 "Propagation curves for aeronautical
 * mobile and radionavigation services using the VHF, UHF and SHF bands", Annex
 * 2, Section 3.
 *
 * @param[in]  d__km       %Path distance, in km
 * @param[in]  h_1__meter  Height of the low terminal, in meters
 * @param[in]  h_2__meter  Height of the high terminal, in meters
 * @param[in]  f__mhz      Frequency, in MHz
 * @param[in]  p           Time percentage
 * @param[in]  T_pol       Polarization type (POLARIZATION__HORIZONTAL or
 *                         POLARIZATION__VERTICAL)
 * @param[out] result      Result structure containing various
 * @param[out] terminal_1  Structure containing parameters dealing with the
 * geometry of terminal 1
 * @param[out] terminal_2  Structure containing parameters dealing with the
 * geometry of terminal 2
 * @param[out] tropo       Structure containing troposcatter parameters
 * @param[out] path        Structure containing propagation path parameters
 * @param[out] los_params  Structure containing line of sight parameters
 * @return                 RETURNCODES value explaining success of code run
 ****************************************************************************/
int P528_Ex(double d__km, double h_1__meter, double h_2__meter, double f__mhz,
            int T_pol, double p, Result* result, Terminal* terminal_1,
            Terminal* terminal_2, TroposcatterParams* tropo, Path* path,
            LineOfSightParams* los_params) {
  // reset Results struct
  result->A_fs__db = 0;
  result->A_a__db = 0;
  result->A__db = 0;
  result->d__km = 0;
  result->theta_h1__rad = 0;
  result->propagation_mode = PROP_MODE__NOT_SET;
  result->warnings = WARNING__NO_WARNINGS;

  int err = ValidateInputs(d__km, h_1__meter, h_2__meter, f__mhz, T_pol, p,
                           &result->warnings);
  if (err != SUCCESS) {
    if (err == ERROR_HEIGHT_AND_DISTANCE) {
      result->A_fs__db = 0;
      result->A_a__db = 0;
      result->A__db = 0;
      result->d__km = 0;
      return SUCCESS;
    } else
      return err;
  }

  // Steps 1-3 and everything else that does not depend on the distance
  ProfileContext context;
  InitializeProfileContext(h_1__meter, h_2__meter, f__mhz, T_pol, p, &context);

  // Steps 4-7 for the requested distance
  EvaluateProfilePoint(&context, d__km, result, tropo, los_params);

  *terminal_1 = context.terminal_1;
  *terminal_2 = context.terminal_2;
  *path = context.path;

  if (result->warnings == WARNING__NO_WARNINGS)
    return SUCCESS;
  else
    return SUCCESS_WITH_WARNINGS;
}

/**
 * Computes everything that does not depend on the path distance: the terminal
 * geometries, the smooth earth diffraction line and the line-of-sight setup.
 * The inputs must have passed ValidateInputs().
 *
 * @param[in]  h_1__meter  Height of the low terminal, in meters
 * @param[in]  h_2__meter  Height of the high terminal, in meters
 * @param[in]  f__mhz      Frequency, in MHz
 * @param[in]  T_pol       Polarization type
 * @param[in]  p           Time percentage
 * @param[out] context     Distance-independent state
 ****************************************************************************/
void InitializeProfileContext(double h_1__meter, double h_2__meter,
                              double f__mhz, int T_pol, double p,
                              ProfileContext* context) {
  Terminal* terminal_1 = &context->terminal_1;
  Terminal* terminal_2 = &context->terminal_2;
  Path* path = &context->path;

  context->f__mhz = f__mhz;
  context->T_pol = T_pol;
  context->p = p;

  /////////////////////////////////////////////
  // Compute terminal geometries
  //

  // Step 1 for low terminal
  terminal_1->h_r__km = h_1__meter / 1000;
  TerminalGeometry(f__mhz, terminal_1);

  // Step 1 for high terminal
  terminal_2->h_r__km = h_2__meter / 1000;
  TerminalGeometry(f__mhz, terminal_2);

  //
  // Compute terminal geometries
  /////////////////////////////////////////////

  // Step 2
  path->d_ML__km = terminal_1->d_r__km + terminal_2->d_r__km;  // [Eqn 3-1]

  /////////////////////////////////////////////
  // Smooth earth diffraction line calculations
  //

  // Step 3.1
  double d_3__km =
      path->d_ML__km + 0.5 * pow(pow(a_e__km, 2) / f__mhz, THIRD);  // [Eqn 3-2]
  double d_4__km =
      path->d_ML__km + 1.5 * pow(pow(a_e__km, 2) / f__mhz, THIRD);  // [Eqn 3-3]

  // Step 3.2
  double A_3__db = SmoothEarthDiffraction(
      terminal_1->d_r__km, terminal_2->d_r__km, f__mhz, d_3__km, T_pol);
  double A_4__db = SmoothEarthDiffraction(
      terminal_1->d_r__km, terminal_2->d_r__km, f__mhz, d_4__km, T_pol);

  // Step 3.3
  double M_d = (A_4__db - A_3__db) / (d_4__km - d_3__km);  // [Eqn 3-4]
  double A_d0 = A_4__db - M_d * d_4__km;                   // [Eqn 3-5]

  // Step 3.4
  double A_dML__db = (M_d * path->d_ML__km) + A_d0;  // [Eqn 3-6]
  path->d_d__km = -(A_d0 / M_d);                     // [Eqn 3-7]

  context->M_d = M_d;
  context->A_d0 = A_d0;
  context->A_dML__db = A_dML__db;

  //
  // End smooth earth diffraction line calculations
  /////////////////////////////////////////////////

  // Line-of-sight quantities shared by all distances. They are also needed
  // beyond the horizon, to get K_LOS.
  LineOfSightSetup(path, terminal_1, terminal_2, f__mhz, -A_dML__db, T_pol,
                   &context->psi_limit, &context->A_d_0__db);
}

/**
 * Computes the distance-independent quantities of the transhorizon region:
 * K_LOS and the crossover between the diffraction and troposcatter models.
 *
 * @param[in,out] context  Distance-independent state
 ****************************************************************************/
static void InitializeTranshorizon(ProfileContext* context) {
  Path* path = &context->path;
  Terminal* terminal_1 = &context->terminal_1;
  Terminal* terminal_2 = &context->terminal_2;

  // get K_LOS
  Result result_K_LOS;
  LineOfSightPoint(path, terminal_1, terminal_2, &context->K_LOS_params,
                   context->f__mhz, -context->A_dML__db, context->p,
                   path->d_ML__km - 1, context->T_pol, context->psi_limit,
                   context->A_d_0__db, &context->los_layers, nullptr,
                   &result_K_LOS, &context->K_LOS);

  // Step 6.  Search past horizon to find crossover point between Diffraction
  // and Troposcatter models
  context->M_d_trans = context->M_d;
  context->A_d0_trans = context->A_d0;
  context->transhorizon_warnings = WARNING__NO_WARNINGS;
  TranshorizonSearch(path, terminal_1, terminal_2, context->f__mhz,
                     context->A_dML__db, &context->M_d_trans,
                     &context->A_d0_trans, &context->d_crx__km, &context->CASE,
                     &context->transhorizon_warnings);

  context->transhorizon_ready = true;
}

/**
 * Computes the basic transmission loss for one distance, reusing the
 * distance-independent state of the profile (Steps 4-7).
 *
 * @param[in,out] context   Distance-independent state
 * @param[in]  d__km        %Path distance, in km
 * @param[in,out] result    Result structure. result->warnings must hold the
 *                          warnings of ValidateInputs() on entry.
 * @param[out] tropo        Structure containing troposcatter parameters
 * @param[out] los_params   Structure containing line of sight parameters
 ****************************************************************************/
void EvaluateProfilePoint(ProfileContext* context, double d__km,
                          Result* result, TroposcatterParams* tropo,
                          LineOfSightParams* los_params) {
  Path* path = &context->path;
  Terminal* terminal_1 = &context->terminal_1;
  Terminal* terminal_2 = &context->terminal_2;
  const double f__mhz = context->f__mhz;
  const double p = context->p;

  // Step 4.  If the path is in the Line-of-Sight range, call LOS and then exit
  if (path->d_ML__km - d__km > 0.001) {
    double K_LOS = 0;
    result->propagation_mode = PROP_MODE__LOS;
    LineOfSightPoint(path, terminal_1, terminal_2, los_params, f__mhz,
                     -context->A_dML__db, p, d__km, context->T_pol,
                     context->psi_limit, context->A_d_0__db,
                     &context->los_layers, context->gamma_table.get(), result,
                     &K_LOS);
    return;
  } else {
    // K_LOS and the crossover search are the same for every distance
    if (!context->transhorizon_ready) InitializeTranshorizon(context);

    *los_params = context->K_LOS_params;
    result->warnings |= context->transhorizon_warnings;

    const double K_LOS = context->K_LOS;
    const double M_d = context->M_d_trans;
    const double A_d0 = context->A_d0_trans;
    const double d_crx__km = context->d_crx__km;
    const int CASE = context->CASE;

    /////////////////////////////////////////////
    // Compute terrain attenuation, A_T__db
    //

    // Step 7.1
    double A_d__db = M_d * d__km + A_d0;  // [Eqn 3-14]

    // Step 7.2
    Troposcatter(terminal_1, terminal_2, d__km, f__mhz, tropo);

    // Step 7.3
    double A_T__db;
    if (d__km < d_crx__km) {
      // always in diffraction if less than d_crx
      A_T__db = A_d__db;
      result->propagation_mode = PROP_MODE__DIFFRACTION;
    } else {
      if (CASE == CASE_1) {
        // select the lower loss mode of propagation
        if (tropo->A_s__db <= A_d__db) {
          A_T__db = tropo->A_s__db;
          result->propagation_mode = PROP_MODE__SCATTERING;
        } else {
          A_T__db = A_d__db;
          result->propagation_mode = PROP_MODE__DIFFRACTION;
        }
      } else  // CASE_2
      {
        A_T__db = tropo->A_s__db;
        result->propagation_mode = PROP_MODE__SCATTERING;
      }
    }

    //
    // Compute terrain attenuation, A_T__db
    /////////////////////////////////////////////

    /////////////////////////////////////////////
    // Compute variability
    //

    // f_theta_h is unity for transhorizon paths
    double f_theta_h = 1;

    // compute the 50% and p% of the long-term variability distribution
    double Y_e__db, Y_e_50__db, dummy;
    LongTermVariability(terminal_1->d_r__km, terminal_2->d_r__km, d__km, f__mhz,
                        p, f_theta_h, -A_T__db, &Y_e__db, &dummy);
    LongTermVariability(terminal_1->d_r__km, terminal_2->d_r__km, d__km, f__mhz,
                        50, f_theta_h, -A_T__db, &Y_e_50__db, &dummy);

    // compute the 50% and p% of the Nakagami-Rice distribution
    double ANGLE = 0.02617993878;  // 1.5 deg
    double K_t__db;
    if (tropo->theta_s >= ANGLE)  // theta_s > 1.5 deg
      K_t__db = 20;
    else if (tropo->theta_s <= 0.0)
      K_t__db = K_LOS;
    else
      K_t__db = (tropo->theta_s * (20.0 - K_LOS) / ANGLE) + K_LOS;

    double Y_pi_50__db = 0.0;  //  zero mean
    double Y_pi__db = NakagamiRice(K_t__db, p);

    // combine the long-term and Nakagami-Rice distributions
    double Y_total__db =
        CombineDistributions(Y_e_50__db, Y_e__db, Y_pi_50__db, Y_pi__db, p);

    //
    // Compute variability
    /////////////////////////////////////////////

    /////////////////////////////////////////////
    // Atmospheric absorption for transhorizon path
    //
    // The ray from the surface to the common volume only depends on h_v, so
    // a profile takes it from a table once that pays off.
    SlantPathAttenuationResult result_v;
    if (!(context->surface_rays &&
          context->surface_rays->Lookup(tropo->h_v__km, &result_v.A_gas__db,
                                        &result_v.a__km)))
      SlantPathAttenuation(f__mhz / 1000, 0, tropo->h_v__km, PI / 2, &result_v,
                           context->gamma_table.get());

    result->A_a__db = terminal_1->A_a__db + terminal_2->A_a__db +
                      2 * result_v.A_gas__db;  // [Eqn 3-17]

    //
    // Atmospheric absorption for transhorizon path
    /////////////////////////////////////////////

    /////////////////////////////////////////////
    // Compute free-space loss
    //

    double r_fs__km = terminal_1->a__km + terminal_2->a__km +
                      2 * result_v.a__km;  // [Eqn 3-18]
    result->A_fs__db =
        20.0 * log10(f__mhz) + 20.0 * log10(r_fs__km) + 32.45;  // [Eqn 3-19]

    //
    // Compute free-space loss
    /////////////////////////////////////////////

    result->d__km = d__km;
    result->A__db = result->A_fs__db + result->A_a__db + A_T__db -
                    Y_total__db;  // [Eqn 3-20]
    result->theta_h1__rad = -terminal_1->theta__rad;
  }
}

}  // namespace p528
}  // namespace PSeries
}  // namespace ITU
}  // namespace ITS

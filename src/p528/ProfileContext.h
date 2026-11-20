/** @file ProfileContext.h
 * Internal header: everything P.528 computes that does not depend on the path
 * distance. It is computed once and shared by all points of a profile.
 */
#pragma once

#include <memory>

#include "p528.h"
#include "p676.h"

namespace ITS {
namespace ITU {
namespace PSeries {
namespace p528 {

/**
 * Distance-independent state of a P.528 calculation, valid for one combination
 * of terminal heights, frequency, polarization and time percentage.
 *
 * A context is not thread safe: parts of it are filled on first use.
 */
struct ProfileContext {
  /** @name Inputs the context was built for */
  ///@{
  double f__mhz;
  int T_pol;
  double p;
  ///@}

  /** @name Terminal geometry and smooth earth diffraction line (Steps 1-3) */
  ///@{
  Terminal terminal_1;
  Terminal terminal_2;
  Path path;
  double M_d;        ///< Slope of the smooth earth diffraction line
  double A_d0;       ///< Intercept of the smooth earth diffraction line
  double A_dML__db;  ///< Diffraction loss at d_ML
  ///@}

  /** @name Line-of-sight region */
  ///@{
  double psi_limit;  ///< Reflection angle separating free space and two-ray
  double A_d_0__db;  ///< Line-of-sight loss at d_0
  /** Atmospheric layers between the two terminals, filled on first use */
  RayTraceLayers los_layers;
  ///@}

  /**
   * Specific attenuation versus height for f__mhz. Optional: if set, rays
   * whose layers change with the distance interpolate the attenuation instead
   * of summing it line by line. That pays off for more than a few points.
   */
  std::unique_ptr<SpecificAttenuationTable> gamma_table;

  /** @name Transhorizon region, filled on first use */
  ///@{
  bool transhorizon_ready = false;
  double K_LOS;                     ///< K-value at the edge of line-of-sight
  LineOfSightParams K_LOS_params;   ///< LOS parameters of the K_LOS evaluation
  double M_d_trans;                 ///< Diffraction line slope after Step 6
  double A_d0_trans;                ///< Diffraction line intercept after Step 6
  double d_crx__km;                 ///< Diffraction/troposcatter crossover
  int CASE;                         ///< Case as defined in Step 6.5
  int transhorizon_warnings;        ///< Warnings raised by the search
  ///@}
};

void InitializeProfileContext(double h_1__meter, double h_2__meter,
                              double f__mhz, int T_pol, double p,
                              ProfileContext* context);

void EvaluateProfilePoint(ProfileContext* context, double d__km,
                          Result* result, TroposcatterParams* tropo,
                          LineOfSightParams* los_params);

void LineOfSightSetup(Path* path, Terminal* terminal_1, Terminal* terminal_2,
                      double f__mhz, double A_dML__db, int T_pol,
                      double* psi_limit, double* A_d_0__db);

void LineOfSightPoint(Path* path, Terminal* terminal_1, Terminal* terminal_2,
                      LineOfSightParams* los_params, double f__mhz,
                      double A_dML__db, double p, double d__km, int T_pol,
                      double psi_limit, double A_d_0__db,
                      RayTraceLayers* los_layers,
                      SpecificAttenuationTable* gamma_table, Result* result,
                      double* K_LOS);

}  // namespace p528
}  // namespace PSeries
}  // namespace ITU
}  // namespace ITS

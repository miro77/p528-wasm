/** @file TerminalGeometry.cpp
 * Computes the terminal geometry
 */
#include <math.h>

#include "p528.h"
#include "p676.h"
#include "p835.h"
#include "ProfileContext.h"

namespace ITS {
namespace ITU {
namespace PSeries {
namespace p528 {

/**
 * Computes the terminal geometry as described
 *
 * References:
 *     - Recommendation ITU‑R P.528‑5 "Propagation curves for aeronautical
 * mobile and radionavigation services using the VHF, UHF and SHF bands", Annex
 * 2, Section 4.
 *
 * @param[in] f__mhz    Frequency, in MHz
 * @param[in] terminal  Structure containing parameters dealing with the
 *                      geometry of the terminal
 ******************************************************************************/
void TerminalGeometry(double f__mhz, Terminal* terminal) {
  TerminalGeometry(f__mhz, terminal, nullptr);
}

/**
 * Computes the terminal geometry, optionally with a table for the specific
 * attenuation of the frequency.
 *
 * With a table only the absorption of the ray to the horizon (A_a__db) is
 * interpolated; it stays within the accuracy of the table. The ray itself is
 * traced with the exact refractive index. So the horizon distance and
 * everything derived from it are bit-identical to those without a table. That
 * matters: they decide where line-of-sight ends and which steps the searches
 * of the line-of-sight region take, and a profile has to agree with P528()
 * on both.
 *
 * @param[in] f__mhz    Frequency, in MHz
 * @param[in] terminal  Structure containing parameters dealing with the
 *                      geometry of the terminal
 * @param[in,out] table Table for f__mhz, or null for the line-by-line
 *                      calculation
 ******************************************************************************/
void TerminalGeometry(double f__mhz, Terminal* terminal,
                      SpecificAttenuationTable* table) {
  double theta_tx__rad = 0;

  SlantPathAttenuationResult result;
  if (table) {
    // SlantPathAttenuation() for an upward ray, with the layers described
    // above
    RayTraceConfig config;
    config.temperature = GlobalTemperature;
    config.dry_pressure = GlobalPressure;
    config.wet_pressure = GlobalWetPressure;

    RayTraceLayers layers;
    GetRayTraceLayers(f__mhz / 1000, 0, terminal->h_r__km, config, &layers,
                      table, true);
    RayTrace(0, PI / 2 - theta_tx__rad, layers, &result);
  } else {
    SlantPathAttenuation(f__mhz / 1000, 0, terminal->h_r__km,
                         PI / 2 - theta_tx__rad, &result);
  }

  terminal->theta__rad = PI / 2 - result.angle__rad;
  terminal->A_a__db = result.A_gas__db;
  terminal->a__km = result.a__km;

  // compute arc distance
  double central_angle = ((PI / 2 - result.angle__rad) - theta_tx__rad +
                          result.bending__rad);  // [Thayer, Equ 2], rearranged
  terminal->d_r__km = a_0__km * central_angle;

  terminal->phi__rad = terminal->d_r__km / a_e__km;  // [Eqn 4-1]

  // [Eqn 4-2]
  terminal->h_e__km = (a_e__km / cos(terminal->phi__rad)) - a_e__km;

  terminal->delta_h__km = terminal->h_r__km - terminal->h_e__km;  // [Eqn 4-3]
}

}  // namespace p528
}  // namespace PSeries
}  // namespace ITU
}  // namespace ITS

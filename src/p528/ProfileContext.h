/** @file ProfileContext.h
 * Internal header: everything P.528 computes that does not depend on the path
 * distance. It is computed once and shared by all points of a profile.
 */
#pragma once

#include <memory>
#include <vector>

#include "p528.h"
#include "p676.h"

namespace ITS {
namespace ITU {
namespace PSeries {
namespace p528 {

/**
 * The reflection angles FindPsiAtDistance() has tried, with their distances.
 *
 * The search halves an interval of angles, so the angles it tries form a
 * binary tree that is the same for every distance. Distances that are close
 * to each other share all but the last steps. The tree stores the upper
 * levels, the lower ones would rarely be reused.
 */
struct PsiSearchTree {
  struct Node {
    double d__km;                ///< Path distance of the angle of this node
    int child[2] = {-1, -1};     ///< Next angle above and below, -1 if not stored
  };
  static constexpr int MAX_DEPTH = 20;     ///< Levels that are stored
  static constexpr size_t MAX_NODES = 1 << 18;  ///< Limits the memory to 4 MB
  std::vector<Node> nodes;                 ///< The first node is the root
};

/**
 * State of the line-of-sight region that all distances share. It is filled on
 * first use.
 */
struct LineOfSightShared {
  /** Atmospheric layers between the two terminals */
  RayTraceLayers layers;
  /** Rays through the layers versus their elevation angle (profiles only) */
  std::unique_ptr<LayerRayTable> rays;
  /** Downward-looking rays between the two terminals (profiles only) */
  std::unique_ptr<GrazingRayCache> grazing_rays;
  /** Reflection angles tried so far (profiles only) */
  PsiSearchTree psi_search;
};

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
  /** Layers and rays between the two terminals, filled on first use */
  LineOfSightShared los;
  ///@}

  /**
   * Specific attenuation versus height for f__mhz. Optional: if set before
   * InitializeProfileContext(), all rays interpolate the attenuation instead
   * of summing it line by line, and only their loss and length are computed
   * where nothing else is needed. That pays off for more than a few points.
   */
  std::unique_ptr<SpecificAttenuationTable> gamma_table;

  /**
   * Rays from the surface up to the common volume, for f__mhz. Optional: if
   * set, the ray of a distance beyond the horizon is interpolated instead of
   * traced. Needs gamma_table.
   */
  std::unique_ptr<SurfaceRayTable> surface_rays;

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

void TerminalGeometry(double f__mhz, Terminal* terminal,
                      SpecificAttenuationTable* table);

double RayOpticsDistance(Terminal* terminal_1, Terminal* terminal_2,
                         double psi);

double FindPsiAtDistance(double d__km, Terminal* terminal_1,
                         Terminal* terminal_2);
double FindPsiAtDistance(double d__km, Terminal* terminal_1,
                         Terminal* terminal_2, PsiSearchTree* tree);

void LineOfSightSetup(Path* path, Terminal* terminal_1, Terminal* terminal_2,
                      double f__mhz, double A_dML__db, int T_pol,
                      double* psi_limit, double* A_d_0__db);

void LineOfSightPoint(Path* path, Terminal* terminal_1, Terminal* terminal_2,
                      LineOfSightParams* los_params, double f__mhz,
                      double A_dML__db, double p, double d__km, int T_pol,
                      double psi_limit, double A_d_0__db,
                      LineOfSightShared* shared,
                      SpecificAttenuationTable* gamma_table, Result* result,
                      double* K_LOS);

}  // namespace p528
}  // namespace PSeries
}  // namespace ITU
}  // namespace ITS

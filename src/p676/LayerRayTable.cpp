#include "p676.h"
#include <cmath>

namespace ITS {
namespace ITU {
namespace PSeries {
namespace p528 {

namespace {

// Positions of the two check points within [-1, 1], away from the nodes
const double CHECK_POINTS[] = { -0.83, 0.37 };

// A segment is built from one ray per node and check point, so it only pays
// off for more rays than that. The first rays are left to the caller.
const int RAYS_BEFORE_BUILD = LayerRayTable::DEGREE + 1 + 2;

// Relative accuracy a segment has to reach to be used
const double TOLERANCE = 1e-10;

// More segments than this would mean layers of a negligible thickness
const int MAX_SEGMENTS = 64;

/*=============================================================================
 |
 |  Description:  Evaluates a Chebyshev series with the Clenshaw recurrence.
 |
 |        Input:  c             - Coefficients c_0 ... c_DEGREE
 |                x             - Position within the segment, in [-1, 1]
 |
 |      Returns:  Value of the series
 |
 *===========================================================================*/
inline double EvaluateChebyshev(const double* c, double x)
{
    double b_1 = 0;
    double b_2 = 0;
    for (int j = LayerRayTable::DEGREE; j >= 1; j--)
    {
        double b_0 = 2 * x * b_1 - b_2 + c[j];
        b_2 = b_1;
        b_1 = b_0;
    }
    return x * b_1 - b_2 + c[0];
}

} // namespace

/*=============================================================================
 |
 |  Description:  Table of the rays through one set of layers, e.g. between
 |                the two terminals of a line-of-sight path. With the layers
 |                fixed, the loss and the length of a ray only depend on
 |                its elevation angle, and smoothly so.
 |
 |                The variable of the table is w = cos(beta_1). Equation 17
 |                makes the path length of layer i a function of
 |                cos^2(beta_i) = k_i + c_i^2 w^2, with constants c_i and
 |                k_i >= 0 of the layer. It is analytic in w for real w,
 |                but k_i is tiny for the lowest layers, so there are
 |                singularities at imaginary w down to about
 |                +-j sqrt(2 delta_1 / r_1) (delta_1: thickness of the
 |                lowest layer). An interpolation in w therefore converges
 |                fast on an interval whose length is comparable to its
 |                distance from w = 0. The table uses one segment from
 |                w = 0 up to half of that bound and one per octave of w
 |                above it. A segment is built when enough rays ask for it.
 |
 |        Input:  h_1__km       - Height at which the rays start, in km
 |                layers        - Layers from GetRayTraceLayers() for
 |                                h_1__km. They have to outlive this object.
 |
 *===========================================================================*/
LayerRayTable::LayerRayTable(double h_1__km, const RayTraceLayers* layers)
    : h_1__km_(h_1__km), layers_(layers)
{
    // distance of the nearest singularity from w = 0
    double delta_1__km = LayerThickness(layers->m, layers->i_lower);
    double w_singular = sqrt(2 * delta_1__km / (a_0__km + h_1__km));

    // No layers, or layers without thickness (both heights equal): there is
    // nothing to interpolate. This also catches NaN.
    if (!(layers->i_upper > layers->i_lower && w_singular > 0 && w_singular < 1))
        return;

    // segment 0 covers [0, 2^exponent_first_), with 2^exponent_first_ <= w_singular / 2
    frexp(w_singular / 2, &exponent_first_);
    exponent_first_ -= 1;
    if (-exponent_first_ + 1 > MAX_SEGMENTS)
        return;

    // segment k > 0 covers the octave [2^(exponent_first_ + k - 1), 2^(exponent_first_ + k)),
    // the last one ends at w = 1
    segments_.resize(-exponent_first_ + 1);

    // discrete cosine transform from node values to Chebyshev coefficients
    const int N = DEGREE + 1;
    for (int k = 0; k < N; k++)
        for (int j = 0; j < N; j++)
            transform_[k][j] = (k == 0 ? 1.0 : 2.0) / N * cos(PI * k * (j + 0.5) / N);
}

/*=============================================================================
 |
 |  Description:  Traces a ray through the layers like RayTraceLoss(), but
 |                for a given cosine of the elevation angle. Unlike the
 |                sine it resolves grazing rays.
 |
 |        Input:  cos_beta_1    - Cosine of the elevation angle from zenith
 |
 |       Output:  A_gas__db     - Gaseous absorption along the ray, in dB
 |                a__km         - Ray length, in km
 |
 *===========================================================================*/
void LayerRayTable::Trace(double cos_beta_1, double* A_gas__db, double* a__km) const
{
    const int i_lower = layers_->i_lower;
    const int i_upper = layers_->i_upper;
    const double m = layers_->m;
    const double exp_lower = LayerExp(i_lower - 1);
    const double exp_step = LayerExp(1) - 1;

    *A_gas__db = 0;
    *a__km = 0;

    double n_i = layers_->n[0];
    double gamma_i = layers_->gamma[0];
    double r_i__km = a_0__km + h_1__km_;
    const double r_1__km = r_i__km;
    const double n_1 = n_i;

    // summation from Equation 13
    for (int i = i_lower; i <= i_upper - 1; i++)
    {
        double h_ii__km = h_1__km_ + m * ((LayerExp(i + 1 - 1) - exp_lower) / exp_step);
        double r_ii__km = a_0__km + h_ii__km;
        double delta_i__km = LayerThickness(m, i);

        // Equation 19b: sin(beta_i) = c_i sin(beta_1), which makes
        // cos^2(beta_i) = (1 - c_i^2) + c_i^2 cos^2(beta_1)
        double c_i = (n_1 * r_1__km) / (n_i * r_i__km);
        double cos2_beta_i = MAX(0, (1 - c_i) * (1 + c_i) + (c_i * c_i) * (cos_beta_1 * cos_beta_1));
        double cos_beta_i = sqrt(cos2_beta_i);

        // path length through ith layer, Equation 17
        double a_i__km = -r_i__km * cos_beta_i + sqrt((r_i__km * r_i__km) * cos2_beta_i + 2 * r_i__km * delta_i__km + delta_i__km * delta_i__km);

        *a__km += a_i__km;
        *A_gas__db += a_i__km * gamma_i;

        // shift for next loop
        n_i = layers_->n[i + 1 - i_lower];
        gamma_i = layers_->gamma[i + 1 - i_lower];
        r_i__km = r_ii__km;
    }
}

/*=============================================================================
 |
 |  Description:  Builds the interpolation of a segment and verifies it. A
 |                segment that does not reach the tolerance is marked as
 |                failed, its rays are traced as before.
 |
 |        Input:  i_segment     - Index of the segment
 |
 *===========================================================================*/
void LayerRayTable::Build(int i_segment)
{
    const int N = DEGREE + 1;
    Segment& segment = segments_[i_segment];

    double lo = (i_segment == 0) ? 0 : ldexp(1.0, exponent_first_ + i_segment - 1);
    double hi = ldexp(1.0, exponent_first_ + i_segment);
    segment.mid = 0.5 * (lo + hi);
    segment.scale = 2 / (hi - lo);

    // sample at the Chebyshev nodes
    double A__db[N];
    double a__km[N];
    double largest__db = 0;
    double largest__km = 0;
    for (int j = 0; j < N; j++)
    {
        double x = cos(PI * (j + 0.5) / N);
        Trace(segment.mid + x / segment.scale, &A__db[j], &a__km[j]);
        largest__db = MAX(largest__db, fabs(A__db[j]));
        largest__km = MAX(largest__km, fabs(a__km[j]));
    }

    for (int k = 0; k < N; k++)
    {
        double sum_A = 0;
        double sum_a = 0;
        for (int j = 0; j < N; j++)
        {
            sum_A += transform_[k][j] * A__db[j];
            sum_a += transform_[k][j] * a__km[j];
        }
        segment.c_A[k] = sum_A;
        segment.c_a[k] = sum_a;
    }

    // A series has converged if its last coefficients are negligible. As an
    // independent check the interpolation is compared with traced rays
    // between the nodes.
    double error__db = fabs(segment.c_A[DEGREE]) + fabs(segment.c_A[DEGREE - 1]);
    double error__km = fabs(segment.c_a[DEGREE]) + fabs(segment.c_a[DEGREE - 1]);
    for (double x : CHECK_POINTS)
    {
        double exact__db, exact__km;
        Trace(segment.mid + x / segment.scale, &exact__db, &exact__km);
        error__db = MAX(error__db, fabs(exact__db - EvaluateChebyshev(segment.c_A, x)));
        error__km = MAX(error__km, fabs(exact__km - EvaluateChebyshev(segment.c_a, x)));
    }

    // both comparisons are false for NaN
    bool accurate = error__db <= TOLERANCE * largest__db && error__km <= TOLERANCE * largest__km;
    segment.state = accurate ? Segment::READY : Segment::FAILED;
}

/*=============================================================================
 |
 |  Description:  Loss and length of the ray that starts with an elevation
 |                angle, like RayTraceLoss().
 |
 |        Input:  beta_1__rad   - Elevation angle (from zenith), in rad
 |
 |       Output:  A_gas__db     - Gaseous absorption along the ray, in dB
 |                a__km         - Ray length, in km
 |
 |      Returns:  true if the outputs are set. Otherwise the ray is not
 |                covered by the table (yet) and has to be traced.
 |
 *===========================================================================*/
bool LayerRayTable::Lookup(double beta_1__rad, double* A_gas__db, double* a__km)
{
    if (segments_.empty())
        return false;

    // The cosine as RayTraceLoss() obtains it for the lowest layer. Near
    // grazing it is limited by the resolution of the sine; using the same
    // value keeps the table consistent with the trace.
    double sin_beta_1 = sin(beta_1__rad);
    double cos_beta_1 = sqrt((1 - sin_beta_1) * (1 + sin_beta_1));

    // upward rays only (this also rejects NaN)
    if (!(beta_1__rad >= 0 && beta_1__rad <= PI / 2 && cos_beta_1 <= 1))
        return false;

    // cos_beta_1 lies in [2^(exponent - 1), 2^exponent)
    int exponent;
    frexp(cos_beta_1, &exponent);
    int i_segment = (cos_beta_1 == 0) ? 0 : exponent - exponent_first_;
    i_segment = MAX(0, MIN(i_segment, (int)segments_.size() - 1));

    Segment& segment = segments_[i_segment];
    if (segment.state != Segment::READY)
    {
        if (segment.state == Segment::FAILED)
            return false;
        if (segment.rays < RAYS_BEFORE_BUILD)
        {
            segment.rays++;
            return false;
        }
        Build(i_segment);
        if (segment.state != Segment::READY)
            return false;
    }

    double x = (cos_beta_1 - segment.mid) * segment.scale;
    *A_gas__db = EvaluateChebyshev(segment.c_A, x);
    *a__km = EvaluateChebyshev(segment.c_a, x);
    return true;
}

}  // namespace p528
}  // namespace PSeries
}  // namespace ITU
}  // namespace ITS

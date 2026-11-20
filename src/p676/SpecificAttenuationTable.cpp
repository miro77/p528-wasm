#include "p676.h"
#include "p835.h"
#include <cmath>

namespace ITS {
namespace ITU {
namespace PSeries {
namespace p528 {

namespace {

// Relative accuracy of the interpolation, for the specific attenuation and
// for n - 1. The line-by-line summation itself is only reproducible to
// roughly 1e-14, so there is nothing to gain below.
const double TOLERANCE = 1e-11;

// Segments are never split below this width. It is only reached at kinks of
// the reference atmosphere, where a narrower segment does not help anymore.
const double MIN_SEGMENT__KM = 1e-6;

// Heights supported by the reference atmosphere of Rec. ITU-R P.835
const double H_MIN__KM = 0;
const double H_MAX__KM = 100;

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
    for (int j = SpecificAttenuationTable::DEGREE; j >= 1; j--)
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
 |  Description:  Table of the specific gaseous attenuation and of the
 |                refractive index versus height, for one frequency.
 |
 |                The line-by-line calculation of Rec. ITU-R P.676 is by far
 |                the most expensive part of a ray trace. For a fixed
 |                frequency it is a function of the height alone, and a
 |                smooth one between the kinks of the reference atmosphere.
 |                The table approximates it by piecewise Chebyshev
 |                interpolation. Segments are created on first use and are
 |                refined until the interpolation error is negligible, so
 |                all rays of a profile share the line-by-line results.
 |                The refractive index is tabulated alongside, it needs
 |                the same atmospheric parameters.
 |
 |        Input:  f__ghz        - Frequency, in GHz
 |                config        - Structure containing atmospheric params
 |
 *===========================================================================*/
SpecificAttenuationTable::SpecificAttenuationTable(double f__ghz, RayTraceConfig config)
    : f__ghz_(f__ghz), config_(config)
{
    // The reference atmosphere is piecewise defined on geopotential height.
    // Start a new cell at each of its boundaries, and at least every km.
    const double regime_boundaries[] = { 11, 20, 32, 47, 51, 71 };
    for (double h_prime__km : regime_boundaries)
        edges_.push_back(ConvertToGeometricHeight(h_prime__km));
    for (int h__km = (int)H_MIN__KM; h__km <= (int)H_MAX__KM; h__km++)
        edges_.push_back(h__km);
    sort(edges_.begin(), edges_.end());

    cells_.resize(edges_.size() - 1);
}

/*=============================================================================
 |
 |  Description:  The exact layer properties at a height.
 |
 |        Input:  h__km         - Height, in km
 |
 |       Output:  n_minus_1     - Refractive index - 1
 |                gamma         - Line-by-line specific gaseous attenuation,
 |                                in dB/km
 |
 *===========================================================================*/
void SpecificAttenuationTable::Exact(double h__km, double* n_minus_1, double* gamma) const
{
    double T__kelvin = config_.temperature(h__km);
    double p__hPa = config_.dry_pressure(h__km);
    double e__hPa = config_.wet_pressure(h__km);

    *n_minus_1 = RefractiveIndexMinusOne(p__hPa, T__kelvin, e__hPa);
    *gamma = SpecificAttenuation(f__ghz_, T__kelvin, e__hPa, p__hPa);
}

/*=============================================================================
 |
 |  Description:  Layer properties at a height, from the table.
 |
 |        Input:  h__km         - Height, in km
 |
 |       Output:  n             - Refractive index
 |                gamma         - Specific gaseous attenuation, in dB/km
 |
 *===========================================================================*/
void SpecificAttenuationTable::Lookup(double h__km, double* n, double* gamma)
{
    // Consecutive lookups are usually neighbouring layers of one ray
    if (!(last_ && h__km >= last_->lo__km && h__km < last_->hi__km))
    {
        // outside of the tabulated range (this includes NaN)
        if (!(h__km >= H_MIN__KM && h__km < H_MAX__KM))
        {
            double n_minus_1;
            Exact(h__km, &n_minus_1, gamma);
            *n = 1 + n_minus_1;
            return;
        }

        last_ = FindSegment(h__km);
    }

    double x = (h__km - last_->mid__km) * last_->scale;
    *n = 1 + EvaluateChebyshev(last_->c_n, x);
    *gamma = EvaluateChebyshev(last_->c_gamma, x);
}

/*=============================================================================
 |
 |  Description:  Finds the segment of a height, creating it if needed.
 |
 |        Input:  h__km         - Height within the tabulated range, in km
 |
 |      Returns:  The segment containing h__km
 |
 *===========================================================================*/
const SpecificAttenuationTable::Segment* SpecificAttenuationTable::FindSegment(double h__km)
{
    // find the cell with edges_[i] <= h__km < edges_[i + 1]
    size_t i_cell = upper_bound(edges_.begin(), edges_.end(), h__km) - edges_.begin() - 1;
    vector<Segment>& cell = cells_[i_cell];
    if (cell.empty())
        Fit(edges_[i_cell], edges_[i_cell + 1], &cell);

    // segments of a cell are sorted by height, find the one containing h__km
    size_t lo = 0;
    size_t hi = cell.size() - 1;
    while (lo < hi)
    {
        size_t mid = (lo + hi) / 2;
        if (h__km < cell[mid].hi__km)
            hi = mid;
        else
            lo = mid + 1;
    }

    return &cell[lo];
}

/*=============================================================================
 |
 |  Description:  Approximates the layer properties between two heights.
 |                The interval is halved recursively until each segment
 |                meets the tolerance. Segments are appended in order of
 |                increasing height.
 |
 |        Input:  lo__km        - Lower height, in km
 |                hi__km        - Upper height, in km
 |
 |       Output:  segments      - Segments covering [lo__km, hi__km)
 |
 *===========================================================================*/
void SpecificAttenuationTable::Fit(double lo__km, double hi__km, vector<Segment>* segments)
{
    const int N = DEGREE + 1;

    Segment segment;
    segment.lo__km = lo__km;
    segment.hi__km = hi__km;
    segment.mid__km = 0.5 * (lo__km + hi__km);
    segment.scale = 2 / (hi__km - lo__km);

    // sample at the Chebyshev nodes
    double gamma[N];
    double n_minus_1[N];
    double largest_gamma = 0;
    double largest_n = 0;
    for (int k = 0; k < N; k++)
    {
        double x = cos(PI * (k + 0.5) / N);
        Exact(segment.mid__km + x / segment.scale, &n_minus_1[k], &gamma[k]);
        largest_gamma = MAX(largest_gamma, fabs(gamma[k]));
        largest_n = MAX(largest_n, fabs(n_minus_1[k]));
    }

    // discrete cosine transform of the samples gives the coefficients
    for (int j = 0; j < N; j++)
    {
        double sum_gamma = 0;
        double sum_n = 0;
        for (int k = 0; k < N; k++)
        {
            double basis = cos(PI * j * (k + 0.5) / N);
            sum_gamma += gamma[k] * basis;
            sum_n += n_minus_1[k] * basis;
        }
        segment.c_gamma[j] = (j == 0 ? 1.0 : 2.0) * sum_gamma / N;
        segment.c_n[j] = (j == 0 ? 1.0 : 2.0) * sum_n / N;
    }

    // A series has converged if its last coefficients are negligible. As an
    // independent check the interpolation is compared with the exact values
    // between the nodes.
    double error_gamma = fabs(segment.c_gamma[DEGREE]) + fabs(segment.c_gamma[DEGREE - 1]);
    double error_n = fabs(segment.c_n[DEGREE]) + fabs(segment.c_n[DEGREE - 1]);
    const double check_points[] = { -0.83, 0.11, 0.67 };
    for (double x : check_points)
    {
        double exact_gamma, exact_n_minus_1;
        Exact(segment.mid__km + x / segment.scale, &exact_n_minus_1, &exact_gamma);
        error_gamma = MAX(error_gamma, fabs(exact_gamma - EvaluateChebyshev(segment.c_gamma, x)));
        error_n = MAX(error_n, fabs(exact_n_minus_1 - EvaluateChebyshev(segment.c_n, x)));
    }

    // both comparisons are false for NaN
    bool accurate = error_gamma <= TOLERANCE * largest_gamma && error_n <= TOLERANCE * largest_n;
    if (accurate || hi__km - lo__km <= MIN_SEGMENT__KM)
    {
        segments->push_back(segment);
    }
    else
    {
        Fit(lo__km, segment.mid__km, segments);
        Fit(segment.mid__km, hi__km, segments);
    }
}

}  // namespace p528
}  // namespace PSeries
}  // namespace ITU
}  // namespace ITS

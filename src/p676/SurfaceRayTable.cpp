#include "p676.h"
#include <cmath>

namespace ITS {
namespace ITU {
namespace PSeries {
namespace p528 {

namespace {

// Rays with fewer layers than this are not tabulated. The layer scale m of
// such short rays varies too much for one interpolation, and tracing them
// directly is cheap.
const int FIRST_LAYERS = 31;

// Layers are added in blocks, as most profiles need the lower ones only
const int BLOCK = 64;

// A ray is traced once per node and check point to build the table, so the
// table only pays off for more rays than that. The first rays are left to
// the caller.
const int RAYS_BEFORE_BUILD = SurfaceRayTable::NODES + 2;

// Positions of the two check points within [-1, 1], away from the nodes
const double CHECK_POINTS[] = { -0.83, 0.37 };

// A layer is interpolated if that reproduces its loss and its length to these
// accuracies. They are well above the rounding noise of a traced layer (about
// 1e-10 dB and 1e-9 km for the lowest layers) and well below the effect of a
// kink of the reference atmosphere (1e-6 dB and more).
const double TOLERANCE__DB = 1e-9;
const double TOLERANCE__KM = 1e-8;

/*=============================================================================
 |
 |  Description:  Evaluates a Chebyshev series with the Clenshaw recurrence.
 |
 |        Input:  c             - Coefficients c_0 ... c_(NODES - 1)
 |                x             - Position within the interval, in [-1, 1]
 |
 |      Returns:  Value of the series
 |
 *===========================================================================*/
inline double EvaluateChebyshev(const double* c, double x)
{
    double b_1 = 0;
    double b_2 = 0;
    for (int j = SurfaceRayTable::NODES - 1; j >= 1; j--)
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
 |  Description:  Table of the rays that leave the surface horizontally and
 |                end at an arbitrary height, for one frequency. P.528 needs
 |                one such ray for every distance beyond the horizon, up to
 |                the height of the common volume.
 |
 |                The loss of such a ray is not a smooth function of its end
 |                height h_2: the number of layers changes with h_2, and with
 |                it the layer scale m of Equation 16c, which shifts every
 |                layer. But for a given m the layers do not depend on h_2 at
 |                all, a higher ray only adds layers on top. And m stays
 |                within a few percent of 1e-4 km, where the contribution of
 |                a layer is a smooth function of m.
 |
 |                The table therefore traces one ray for each of a few
 |                values of m (Chebyshev nodes) and stores, layer by layer,
 |                the interpolation coefficients of the sums up to that
 |                layer. A ray is then a single interpolation in m.
 |
 |                Where the midpoint of a layer crosses a kink of the
 |                reference atmosphere within the range of m, its
 |                contribution is not smooth. Such layers are detected by
 |                their interpolation error, left out of the sums and
 |                evaluated for each ray instead.
 |
 |        Input:  table         - Layer properties for the frequency. It
 |                                has to outlive this object.
 |
 *===========================================================================*/
SurfaceRayTable::SurfaceRayTable(SpecificAttenuationTable* table)
    : table_(table)
{
    // Range of m for rays with more than FIRST_LAYERS layers. The largest
    // value belongs to a ray that ends exactly on a layer boundary of the
    // grid, the smallest to one that ends just above it (Equation 16).
    const int i_upper_min = FIRST_LAYERS + 1;
    double m_max = 1e-4;
    double m_min = 1e-4 * (LayerExp(i_upper_min - 2) - 1) / (LayerExp(i_upper_min - 1) - 1);
    m_mid_ = 0.5 * (m_min + m_max);
    m_half_ = 0.5 * (m_max - m_min);

    for (int j = 0; j < NODES; j++)
        m_node_[j] = m_mid_ + m_half_ * cos(PI * (j + 0.5) / NODES);
    for (int j = 0; j < CHECKS; j++)
        m_node_[NODES + j] = m_mid_ + m_half_ * CHECK_POINTS[j];

    // discrete cosine transform from node values to Chebyshev coefficients
    for (int k = 0; k < NODES; k++)
        for (int j = 0; j < NODES; j++)
            transform_[k][j] = (k == 0 ? 1.0 : 2.0) / NODES * cos(PI * k * (j + 0.5) / NODES);

    // no layers yet: the sums are zero
    c_A_.assign(NODES, 0);
    c_a_.assign(NODES, 0);
}

/*=============================================================================
 |
 |  Description:  Path length and loss of one layer of a ray that leaves the
 |                surface horizontally, as in GetRayTraceLayers() and
 |                RayTraceLoss() for h_1 = 0 km and beta_1 = 90 degrees.
 |
 |        Input:  i             - Layer, the lowest one is layer 1
 |                m             - Layer scale of the ray, Equation 16c
 |                n_1           - Refractive index of layer 1
 |
 |       Output:  a_i__km       - Path length through the layer, in km
 |                A_i__db       - Loss of the layer, in dB
 |
 *===========================================================================*/
void SurfaceRayTable::LayerTerms(int i, double m, double n_1, double* a_i__km, double* A_i__db) const
{
    const double sin_beta_1 = sin(PI / 2);
    const double r_1__km = a_0__km;

    double h_i__km = m * ((LayerExp(i - 1) - LayerExp(0)) / (LayerExp(1) - 1));
    double delta_i__km = LayerThickness(m, i);
    double r_i__km = a_0__km + h_i__km;

    double n_i, gamma_i;
    table_->Lookup(h_i__km + delta_i__km / 2, &n_i, &gamma_i);

    // Equations 19b and 17
    double sin_beta_i = MIN(1, (n_1 * r_1__km) / (n_i * r_i__km) * sin_beta_1);
    double cos_beta_i = sqrt((1 - sin_beta_i) * (1 + sin_beta_i));
    *a_i__km = -r_i__km * cos_beta_i + sqrt(pow(r_i__km, 2) * pow(cos_beta_i, 2) + 2 * r_i__km * delta_i__km + pow(delta_i__km, 2));
    *A_i__db = *a_i__km * gamma_i;
}

/*=============================================================================
 |
 |  Description:  Extends the table to more layers.
 |
 |        Input:  n_layers      - Number of layers the table has to cover
 |
 *===========================================================================*/
void SurfaceRayTable::Extend(int n_layers)
{
    if (n_layers_ == 0)
    {
        double gamma;
        for (int j = 0; j < NODES + CHECKS; j++)
            table_->Lookup(m_node_[j] / 2, &n_1_[j], &gamma);
    }

    c_A_.resize((size_t)(n_layers + 1) * NODES);
    c_a_.resize((size_t)(n_layers + 1) * NODES);

    for (int i = n_layers_ + 1; i <= n_layers; i++)
    {
        double a__km[NODES + CHECKS];
        double A__db[NODES + CHECKS];
        for (int j = 0; j < NODES + CHECKS; j++)
            LayerTerms(i, m_node_[j], n_1_[j], &a__km[j], &A__db[j]);

        double c_A[NODES];
        double c_a[NODES];
        for (int k = 0; k < NODES; k++)
        {
            double sum_A = 0;
            double sum_a = 0;
            for (int j = 0; j < NODES; j++)
            {
                sum_A += transform_[k][j] * A__db[j];
                sum_a += transform_[k][j] * a__km[j];
            }
            c_A[k] = sum_A;
            c_a[k] = sum_a;
        }

        // The interpolation has converged if its last coefficients are
        // negligible. As an independent check it is compared with the layer
        // itself between the nodes.
        double error__db = fabs(c_A[NODES - 1]) + fabs(c_A[NODES - 2]);
        double error__km = fabs(c_a[NODES - 1]) + fabs(c_a[NODES - 2]);
        for (int j = 0; j < CHECKS; j++)
        {
            error__db = MAX(error__db, fabs(EvaluateChebyshev(c_A, CHECK_POINTS[j]) - A__db[NODES + j]));
            error__km = MAX(error__km, fabs(EvaluateChebyshev(c_a, CHECK_POINTS[j]) - a__km[NODES + j]));
        }

        // both comparisons are false for NaN
        bool smooth = error__db <= TOLERANCE__DB && error__km <= TOLERANCE__KM;
        if (!smooth)
            irregular_.push_back(i);

        const double* previous_A = &c_A_[(size_t)(i - 1) * NODES];
        const double* previous_a = &c_a_[(size_t)(i - 1) * NODES];
        double* sums_A = &c_A_[(size_t)i * NODES];
        double* sums_a = &c_a_[(size_t)i * NODES];
        for (int k = 0; k < NODES; k++)
        {
            sums_A[k] = previous_A[k] + (smooth ? c_A[k] : 0);
            sums_a[k] = previous_a[k] + (smooth ? c_a[k] : 0);
        }
    }

    n_layers_ = n_layers;
}

/*=============================================================================
 |
 |  Description:  Loss and length of the ray from the surface (h_1 = 0 km,
 |                launched horizontally) up to a height.
 |
 |        Input:  h_2__km       - Height at which the ray ends, in km
 |
 |       Output:  A_gas__db     - Gaseous absorption along the ray, in dB
 |                a__km         - Ray length, in km
 |
 |      Returns:  true if the outputs are set. Otherwise the ray is not
 |                covered by the table (yet) and has to be traced.
 |
 *===========================================================================*/
bool SurfaceRayTable::Lookup(double h_2__km, double* A_gas__db, double* a__km)
{
    // the reference atmosphere ends at 100 km (this also rejects NaN)
    if (!(h_2__km > 0 && h_2__km < 100))
        return false;

    // Equations 16(a)-(c), as in GetRayTraceLayers() for h_1 = 0 km
    const int i_lower = 1;
    int i_upper = (int)ceil(100 * log(1e4 * h_2__km * (exp(1. / 100.) - 1) + 1) + 1);
    int n_layers = i_upper - i_lower;
    if (n_layers < FIRST_LAYERS || n_layers > LAST_LAYERS)
        return false;

    double m = ((LayerExp(2) - LayerExp(1)) / (LayerExp(i_upper) - LayerExp(i_lower))) * h_2__km;

    // rounding may put m slightly outside of its range, more is not expected
    double x = (m - m_mid_) / m_half_;
    if (!(fabs(x) <= 1 + 1e-6))
        return false;

    if (rays_ < RAYS_BEFORE_BUILD)
    {
        rays_++;
        return false;
    }
    if (n_layers > n_layers_)
        Extend(MIN(LAST_LAYERS, n_layers + BLOCK));

    *A_gas__db = EvaluateChebyshev(&c_A_[(size_t)n_layers * NODES], x);
    *a__km = EvaluateChebyshev(&c_a_[(size_t)n_layers * NODES], x);

    // layers that are not part of the sums
    if (!irregular_.empty() && irregular_.front() <= n_layers)
    {
        double n_1, gamma_1;
        table_->Lookup(m / 2, &n_1, &gamma_1);
        for (int i : irregular_)
        {
            if (i > n_layers)
                break;
            double a_i__km, A_i__db;
            LayerTerms(i, m, n_1, &a_i__km, &A_i__db);
            *A_gas__db += A_i__db;
            *a__km += a_i__km;
        }
    }

    return true;
}

}  // namespace p528
}  // namespace PSeries
}  // namespace ITU
}  // namespace ITS

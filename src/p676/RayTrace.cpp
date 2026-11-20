#include "p676.h"
#include "p835.h"
#include <cmath>

namespace ITS {
namespace ITU {
namespace PSeries {
namespace p528 {

/*=============================================================================
 |
 |  Description:  exp(k / 100), the growth of the layer thickness. A ray
 |                trace needs it several times per layer, always for a
 |                small integer k, so the values are tabulated once. The
 |                table holds the results of the very same expression,
 |                which keeps the ray trace bit-identical.
 |
 |        Input:  k             - Layer index
 |
 |      Returns:  exp(k / 100)
 |
 *===========================================================================*/
double LayerExp(int k)
{
    // the top of the reference atmosphere (100 km) is at layer 922
    static const int SIZE = 1024;
    static const vector<double> table = [] {
        vector<double> values(SIZE);
        for (int j = 0; j < SIZE; j++)
            values[j] = exp(j / 100.);
        return values;
    }();

    if (k >= 0 && k < SIZE)
        return table[k];
    return exp(k / 100.);
}

/*=============================================================================
 |
 |  Description:  Thickness of the ith layer.
 |
 |        Input:  m             - Internal parameter
 |                i             - Layer of interest
 |
 |      Returns:  delta_i__km   - Layer thickness, in km
 |
 *===========================================================================*/
double LayerThickness(double m, int i)
{
    // Equation 14
    double delta_i__km = m * LayerExp(i - 1);

    return delta_i__km;
}

/*=============================================================================
 |
 |  Description:  Determines the layers between h_1 and h_2 and their
 |                properties (refractive index and specific attenuation).
 |                These depend on the frequency and the two heights only,
 |                not on the elevation angle. Rays that share both heights
 |                (e.g. the points of a profile) can therefore reuse them.
 |
 |        Input:  f__ghz        - Frequency, in GHz
 |                h_1__km       - Height of the low terminal, in km
 |                h_2__km       - Height of the high terminal, in km
 |                config        - Structure containing atmospheric params
 |                table         - Optional table for the frequency. If given,
 |                                the layer properties are interpolated
 |                                instead of computed for each layer.
 |
 |       Output:  layers        - Layer structure
 |
 |      Returns:  [void]
 |
 *===========================================================================*/
static void GetLayerProperties(double f__ghz, double h_i__km, RayTraceConfig config,
    SpecificAttenuationTable* table, double* n, double* gamma)
{
    if (!table)
    {
        GetLayerProperties(f__ghz, h_i__km, config, n, gamma);
        return;
    }

    table->Lookup(h_i__km, n, gamma);
}

void GetRayTraceLayers(double f__ghz, double h_1__km, double h_2__km,
    RayTraceConfig config, RayTraceLayers* layers, SpecificAttenuationTable* table)
{
    // Equations 16(a)-(c)
    int i_lower = floor(100 * log(1e4 * h_1__km * (exp(1. / 100.) - 1) + 1) + 1);
    int i_upper = ceil(100 * log(1e4 * h_2__km * (exp(1. / 100.) - 1) + 1) + 1);
    double m = ((exp(2. / 100.) - exp(1. / 100.)) / (exp(i_upper / 100.) - exp(i_lower / 100.))) * (h_2__km - h_1__km);

    layers->i_lower = i_lower;
    layers->i_upper = i_upper;
    layers->m = m;

    // layers i_lower ... i_upper are stored at index 0 ... i_upper - i_lower
    int count = MAX(i_upper - i_lower + 1, 1);
    layers->n.resize(count);
    layers->gamma.resize(count);

    // starting layer
    double delta_i__km = LayerThickness(m, i_lower);
    double h_i__km = h_1__km + m * ((LayerExp(i_lower - 1) - LayerExp(i_lower - 1)) / (LayerExp(1) - 1));
    GetLayerProperties(f__ghz, h_i__km + delta_i__km / 2, config, table, &layers->n[0], &layers->gamma[0]);

    for (int i = i_lower; i <= i_upper - 1; i++)
    {
        double delta_ii__km = LayerThickness(m, i + 1);
        double h_ii__km = h_1__km + m * ((LayerExp(i + 1 - 1) - LayerExp(i_lower - 1)) / (LayerExp(1) - 1));

        GetLayerProperties(f__ghz, h_ii__km + delta_ii__km / 2, config, table,
            &layers->n[i + 1 - i_lower], &layers->gamma[i + 1 - i_lower]);
    }

    layers->ready = true;
}

/*=============================================================================
 |
 |  Description:  Traces the ray from terminal h_1 to terminal h_2 and
 |                computes results such as atmospheric absorption loss and
 |                ray path length.
 |
 |        Input:  f__ghz        - Frequency, in GHz
 |                h_1__km       - Height of the low terminal, in km
 |                h_2__km       - Height of the high terminal, in km
 |                beta_1__rad   - Elevation angle (from zenith), in rad
 |                config        - Structure containing atmospheric params
 |                table         - Optional table for the frequency. If given,
 |                                the fast path for profiles is used: layer
 |                                properties are interpolated and only
 |                                A_gas__db, a__km and delta_L__km are
 |                                computed, see RayTraceLoss().
 |
 |       Output:  result        - Ray trace result structure
 |
 |      Returns:  [void]
 |
 *===========================================================================*/
void RayTrace(double f__ghz, double h_1__km, double h_2__km, double beta_1__rad,
    RayTraceConfig config, SlantPathAttenuationResult* result,
    SpecificAttenuationTable* table)
{
    RayTraceLayers layers;
    GetRayTraceLayers(f__ghz, h_1__km, h_2__km, config, &layers, table);
    if (table)
        RayTraceLoss(h_1__km, beta_1__rad, layers, result);
    else
        RayTrace(h_1__km, beta_1__rad, layers, result);
}

/*=============================================================================
 |
 |  Description:  Traces the ray through precomputed layers, like RayTrace(),
 |                but only for the path length and the losses. Without the
 |                bending and the arrival angle no inverse trigonometric
 |                function is needed per layer, which makes the trace
 |                several times faster. Compared with RayTrace() the
 |                results differ by rounding only.
 |
 |        Input:  h_1__km       - Height of the low terminal, in km
 |                beta_1__rad   - Elevation angle (from zenith), in rad
 |                layers        - Layers from GetRayTraceLayers() for the
 |                                same h_1__km
 |
 |       Output:  result        - Ray trace result structure. A_gas__db,
 |                                a__km and delta_L__km are set,
 |                                bending__rad and angle__rad are NaN.
 |
 |      Returns:  [void]
 |
 *===========================================================================*/
void RayTraceLoss(double h_1__km, double beta_1__rad, const RayTraceLayers& layers,
    SlantPathAttenuationResult* result)
{
    const int i_lower = layers.i_lower;
    const int i_upper = layers.i_upper;
    const double m = layers.m;

    result->A_gas__db = 0;
    result->a__km = 0;
    result->delta_L__km = 0;
    result->bending__rad = NAN;
    result->angle__rad = NAN;

    // starting layer
    double h_i__km = h_1__km + m * ((LayerExp(i_lower - 1) - LayerExp(i_lower - 1)) / (LayerExp(1) - 1));
    double n_i = layers.n[0];
    double gamma_i = layers.gamma[0];
    double r_i__km = a_0__km + h_i__km;

    // bottom layer properties and launch angle, constant along the ray
    const double r_1__km = r_i__km;
    const double n_1 = n_i;
    const double sin_beta_1 = sin(beta_1__rad);

    // summation from Equation 13
    for (int i = i_lower; i <= i_upper - 1; i++)
    {
        double h_ii__km = h_1__km + m * ((LayerExp(i + 1 - 1) - LayerExp(i_lower - 1)) / (LayerExp(1) - 1));
        double r_ii__km = a_0__km + h_ii__km;
        double delta_i__km = LayerThickness(m, i);

        // Equation 19b gives sin(beta_i). The path length only needs
        // cos(beta_i); the product form stays accurate for grazing rays.
        double sin_beta_i = MIN(1, (n_1 * r_1__km) / (n_i * r_i__km) * sin_beta_1);
        double cos_beta_i = sqrt((1 - sin_beta_i) * (1 + sin_beta_i));

        // path length through ith layer, Equation 17
        double a_i__km = -r_i__km * cos_beta_i + sqrt(pow(r_i__km, 2) * pow(cos_beta_i, 2) + 2 * r_i__km * delta_i__km + pow(delta_i__km, 2));

        result->a__km += a_i__km;
        result->A_gas__db += a_i__km * gamma_i;
        result->delta_L__km += a_i__km * (n_i - 1);     // summation, Equation 23

        // shift for next loop
        n_i = layers.n[i + 1 - i_lower];
        gamma_i = layers.gamma[i + 1 - i_lower];
        r_i__km = r_ii__km;
    }
}

/*=============================================================================
 |
 |  Description:  Traces the ray through precomputed layers. This is the
 |                angle-dependent (and cheap) part of the ray trace.
 |
 |        Input:  h_1__km       - Height of the low terminal, in km
 |                beta_1__rad   - Elevation angle (from zenith), in rad
 |                layers        - Layers from GetRayTraceLayers() for the
 |                                same h_1__km
 |
 |       Output:  result        - Ray trace result structure
 |
 |      Returns:  [void]
 |
 *===========================================================================*/
void RayTrace(double h_1__km, double beta_1__rad, const RayTraceLayers& layers,
    SlantPathAttenuationResult* result)
{
    const int i_lower = layers.i_lower;
    const int i_upper = layers.i_upper;
    const double m = layers.m;

    double gamma_i;
    double gamma_ii;
    double n_i;
    double n_ii;
    double r_i__km;
    double r_ii__km;
    double h_i__km;
    double h_ii__km;
    double a_i__km;
    double alpha_i__rad = beta_1__rad;
    double delta_i__km;
    double beta_i__rad;
    double beta_ii__rad = beta_1__rad;

    // initialize results
    result->A_gas__db = 0;
    result->bending__rad = 0;
    result->a__km = 0;
    result->delta_L__km = 0;

    // initialize starting layer
    h_i__km = h_1__km + m * ((LayerExp(i_lower - 1) - LayerExp(i_lower - 1)) / (LayerExp(1) - 1));
    n_i = layers.n[0];
    gamma_i = layers.gamma[0];
    r_i__km = a_0__km + h_i__km;

    // record bottom layer properties for alpha and beta calculations
    double r_1__km = r_i__km;
    double n_1 = n_i;

    const double sin_beta_1 = sin(beta_1__rad);

    // summation from Equation 13
    for (int i = i_lower; i <= i_upper - 1; i++)
    {
        h_ii__km = h_1__km + m * ((LayerExp(i + 1 - 1) - LayerExp(i_lower - 1)) / (LayerExp(1) - 1));

        n_ii = layers.n[i + 1 - i_lower];
        gamma_ii = layers.gamma[i + 1 - i_lower];

        r_ii__km = a_0__km + h_ii__km;

        delta_i__km = LayerThickness(m, i);

        // Equation 19b
        beta_i__rad = asin(MIN(1, (n_1 * r_1__km) / (n_i * r_i__km) * sin_beta_1));

        // entry angle into the layer interface, Equation 18a
        alpha_i__rad = asin(MIN(1, (n_1 * r_1__km) / (n_i * r_ii__km) * sin_beta_1));

        // path length through ith layer, Equation 17
        const double cos_beta_i = cos(beta_i__rad);
        a_i__km = -r_i__km * cos_beta_i + sqrt(pow(r_i__km, 2) * pow(cos_beta_i, 2) + 2 * r_i__km * delta_i__km + pow(delta_i__km, 2));

        result->a__km += a_i__km;
        result->A_gas__db += a_i__km * gamma_i;
        result->delta_L__km += a_i__km * (n_i - 1);     // summation, Equation 23

        beta_ii__rad = asin(n_i / n_ii * sin(alpha_i__rad));

        // summation of the bending angle, Equation 22a
        // the summation only goes to i_max - 1
        if (i != i_upper - 1)
            result->bending__rad += beta_ii__rad - alpha_i__rad;

        h_i__km = h_ii__km;
        n_i = n_ii;
        gamma_i = gamma_ii;
        r_i__km = r_ii__km;
    }

    result->angle__rad = alpha_i__rad;
}

/*=============================================================================
 |
 |  Description:  Determine the parameters for the ith layer
 |
 |        Input:  f__ghz        - Frequency, in GHz
 |                h_i__km       - Height of the ith layer, in km
 |                config        - Structure containing atmospheric params
 |
 |       Output:  n             - Refractive index
 |                gamma         - Specific attenuation, in dB/km
 |
 |      Returns:  [void]
 |
 *===========================================================================*/
void GetLayerProperties(double f__ghz, double h_i__km, RayTraceConfig config,
    double* n, double* gamma)
{
    // use function pointers to get atmospheric parameters
    double T__kelvin = config.temperature(h_i__km);
    double p__hPa = config.dry_pressure(h_i__km);
    double e__hPa = config.wet_pressure(h_i__km);

    // compute the refractive index for the current layer
    *n = RefractiveIndex(p__hPa, T__kelvin, e__hPa);

    // specific attenuation of layer
    *gamma = SpecificAttenuation(f__ghz, T__kelvin, e__hPa, p__hPa);
}
}  // namespace p528
}  // namespace PSeries
}  // namespace ITU
}  // namespace ITS

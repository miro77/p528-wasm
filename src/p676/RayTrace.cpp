#include "p676.h"
#include "p835.h"

#include <math.h>
#include <numeric>

#include <experimental/simd>
namespace stdx = std::experimental;

using doublev = stdx::native_simd<double>;

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
    double delta_i__km = m * exp((i - 1) / 100.);

    return delta_i__km;
}

doublev LayerThicknessSimd(double m, doublev i)
{
    // Equation 14
    doublev delta_i__km = m * exp((i - 1) / 100.);

    return delta_i__km;
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
 |
 |       Output:  result        - Ray trace result structure
 |
 |      Returns:  [void]
 |
 *===========================================================================*/
void RayTrace(double f__ghz, double h_1__km, double h_2__km, double beta_1__rad,
              RayTraceConfig config, SlantPathAttenuationResult *result)
{
    // Equations 16(a)-(c)
    int i_lower = floor(100 * log(1e4 * h_1__km * (exp(1. / 100.) - 1) + 1) + 1);
    int i_upper = ceil(100 * log(1e4 * h_2__km * (exp(1. / 100.) - 1) + 1) + 1);
    double m = ((exp(2. / 100.) - exp(1. / 100.)) / (exp(i_upper / 100.) - exp(i_lower / 100.))) * (h_2__km - h_1__km);

    // initialize results
    result->A_gas__db = 0.0;
    result->bending__rad = 0.0;
    result->a__km = 0.0;
    result->delta_L__km = 0;
    result->angle__rad = beta_1__rad;

    // initialize starting layer
    double delta_i_lower_km = LayerThickness(m, i_lower);
    double h_i__km = h_1__km + m * ((exp((i_lower - 1) / 100.) - exp((i_lower - 1) / 100.)) / (exp(1 / 100.) - 1));
    double gamma_i;
    double n_i;
    GetLayerProperties(f__ghz, h_i__km + delta_i_lower_km / 2, config, &n_i, &gamma_i);
    double r_i__km = a_0__km + h_i__km;

    // record bottom layer properties for alpha and beta calculations
    double r_1__km = r_i__km;
    double n_1 = n_i;

    const size_t loop_size = i_upper - i_lower + 1; // Größe des Vektors
    std::vector<double> delta_ii__km_vector{delta_i_lower_km};
    std::vector<double> h_ii__km_vector{h_i__km};
    // summation from Equation 13
    delta_ii__km_vector.reserve(loop_size);
    h_ii__km_vector.reserve(loop_size);

    std::vector<int> loop_index_vector(loop_size);
    std::iota(loop_index_vector.begin(), loop_index_vector.end(), i_lower);

    // summation from Equation 13
    int vector_index{0};
    std::size_t i_simd = i_lower;
    for (; i_simd + doublev::size() <= i_upper; i_simd += doublev::size())
    {
        const doublev loop_index_simd = doublev(&loop_index_vector[vector_index + 1], stdx::element_aligned);

        doublev delta_ii__km = LayerThicknessSimd(m, loop_index_simd);
        for (size_t i = 0; i < doublev::size(); i++)
        {
            delta_ii__km_vector.push_back(delta_ii__km[i]);
        }

        doublev h_ii__km = h_1__km + m * ((exp((loop_index_simd - 1) / 100.) - exp((i_lower - 1) / 100.)) / (exp(1 / 100.) - 1));
        for (size_t i = 0; i < doublev::size(); i++)
        {
            h_ii__km_vector.push_back(h_ii__km[i]);
        }
        vector_index += doublev::size();
    }
    for (; i_simd <= i_upper - 1; i_simd++)
    {
        double delta_ii__km = LayerThickness(m, i_simd + 1);
        delta_ii__km_vector.push_back(delta_ii__km);
        double h_ii__km = h_1__km + m * ((exp((i_simd + 1 - 1) / 100.) - exp((i_lower - 1) / 100.)) / (exp(1 / 100.) - 1));
        h_ii__km_vector.push_back(h_ii__km);
    }

    std::vector<double> n_ii_vector{n_i};
    std::vector<double> gamma_ii_vector{gamma_i};
    n_ii_vector.reserve(loop_size);
    gamma_ii_vector.reserve(loop_size);
    vector_index = 0;
    for (int i = i_lower; i <= i_upper - 1; i++)
    {

        double n_ii, gamma_ii;
        GetLayerProperties(f__ghz, h_ii__km_vector[vector_index + 1] + delta_ii__km_vector[vector_index + 1] / 2, config, &n_ii, &gamma_ii);
        n_ii_vector.push_back(n_ii);
        gamma_ii_vector.push_back(gamma_ii);
        ++vector_index;
    }

    vector_index = 0;
    i_simd = i_lower;
    for (; i_simd + doublev::size() <= i_upper; i_simd += doublev::size())
    {
        const doublev h_i__km_simd = doublev(&h_ii__km_vector[vector_index], stdx::element_aligned);
        doublev r_i__km = a_0__km + h_i__km_simd;

        const doublev h_ii__km_simd = doublev(&h_ii__km_vector[vector_index + 1], stdx::element_aligned);
        doublev r_ii__km = a_0__km + h_ii__km_simd;

        doublev delta_i__km = doublev(&delta_ii__km_vector[vector_index], stdx::element_aligned);

        const doublev n_i_simd = doublev(&n_ii_vector[vector_index], stdx::element_aligned);
        const doublev n_ii_simd = doublev(&n_ii_vector[vector_index + 1], stdx::element_aligned);

        // Equation 19b
        doublev beta_i__rad = asin(stdx::min(doublev{1}, (n_1 * r_1__km) / (n_i_simd * r_i__km) * sin(beta_1__rad)));
        // entry angle into the layer interface, Equation 18a
        doublev alpha_i__rad = asin(stdx::min(doublev{1}, (n_1 * r_1__km) / (n_i_simd * r_ii__km) * sin(beta_1__rad)));

        // path length through ith layer, Equation 17
        doublev a_i__km = -r_i__km * cos(beta_i__rad) + sqrt(pow(r_i__km, 2) * pow(cos(beta_i__rad), 2) + 2 * r_i__km * delta_i__km + pow(delta_i__km, 2));

        result->a__km += stdx::reduce(a_i__km);
        const doublev gamma_i_simd = doublev(&gamma_ii_vector[vector_index], stdx::element_aligned);
        result->A_gas__db += stdx::reduce(a_i__km * gamma_i_simd);
        result->delta_L__km += stdx::reduce(a_i__km * (n_i_simd - 1)); // summation, Equation 23

        doublev beta_ii__rad = asin(n_i_simd / n_ii_simd * sin(alpha_i__rad));

        // summation of the bending angle, Equation 22a
        // the summation only goes to i_max - 1
        for (size_t i = 0; i < doublev::size(); i++)
        {
            result->angle__rad = alpha_i__rad[i];
            if (i_simd + i != i_upper - 1)
                result->bending__rad += beta_ii__rad[i] - alpha_i__rad[i];
        }

        vector_index += doublev::size();
    }

    for (; i_simd <= i_upper - 1; i_simd++)
    {
        double r_i__km = a_0__km + h_ii__km_vector[vector_index];
        double r_ii__km = a_0__km + h_ii__km_vector[vector_index + 1];

        double delta_i__km = delta_ii__km_vector[vector_index];

        // Equation 19b
        double beta_i__rad = asin(MIN(1, (n_1 * r_1__km) / (n_ii_vector[vector_index] * r_i__km) * sin(beta_1__rad)));
        // entry angle into the layer interface, Equation 18a
        double alpha_i__rad = asin(MIN(1, (n_1 * r_1__km) / (n_ii_vector[vector_index] * r_ii__km) * sin(beta_1__rad)));

        // path length through ith layer, Equation 17
        double a_i__km = -r_i__km * cos(beta_i__rad) + sqrt(pow(r_i__km, 2) * pow(cos(beta_i__rad), 2) + 2 * r_i__km * delta_i__km + pow(delta_i__km, 2));

        result->a__km += a_i__km;
        result->A_gas__db += a_i__km * gamma_ii_vector[vector_index];
        result->delta_L__km += a_i__km * (n_ii_vector[vector_index] - 1); // summation, Equation 23

        double beta_ii__rad = asin(n_ii_vector[vector_index] / n_ii_vector[vector_index + 1] * sin(alpha_i__rad));

        result->angle__rad = alpha_i__rad;
        // summation of the bending angle, Equation 22a
        // the summation only goes to i_max - 1
        if (i_simd != i_upper - 1)
            result->bending__rad += beta_ii__rad - alpha_i__rad;

        // // shift for next loop
        // h_i__km = h_ii__km_vector[vector_index];
        // n_i = n_ii_vector[vector_index];
        // gamma_i = gamma_ii_vector[vector_index];
        // r_i__km = r_ii__km;

        ++vector_index;
    }
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
                        double *n, double *gamma)
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
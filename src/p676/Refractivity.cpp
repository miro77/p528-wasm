#include "p676.h"

#include <math.h>

#include <experimental/simd>
namespace stdx = std::experimental;

using doublev = stdx::native_simd<double>;

/*=============================================================================
 |
 |  Description:  Imaginary part of the frequency-dependent complex
 |                refractivity due to oxygen.  See Equation (2a).
 |
 |        Input:  f__ghz        - Frequency, in GHz
 |                T__kelvin     - Temperature, in Kelvin
 |                e__hPa        - Water vapour partial pressure, in hPa
 |                p__hPa        - Dry air pressure, in hPa
 |
 |      Returns:  N_o           - Refractivity, in N-Units
 |
 *===========================================================================*/
double OxygenRefractivity(double f__ghz, double T__kelvin, double e__hPa, double p__hPa)
{
    double theta = 300 / T__kelvin;

    double N = 0;

    for (int i = 0; i < OxygenData::f_0.size(); i++)
    {
        // Equation 3, for oxygen
        double S_i = OxygenData::a_1[i] * 1e-7 * p__hPa * pow(theta, 3) * exp(OxygenData::a_2[i] * (1 - theta));

        // compute the width of the line, Equation 6a, for oxygen
        double delta_f__ghz = OxygenData::a_3[i] * 1e-4 * (p__hPa * pow(theta, (0.8 - OxygenData::a_4[i])) + 1.1 * e__hPa * theta);

        // modify the line width to account for Zeeman splitting of the oxygen lines
        // Equation 6b, for oxygen
        delta_f__ghz = sqrt(pow(delta_f__ghz, 2) + 2.25e-6);

        // correction factor due to interference effects in oxygen lines
        // Equation 7, for oxygen
        double delta = (OxygenData::a_5[i] + OxygenData::a_6[i] * theta) * 1e-4 * (p__hPa + e__hPa) * pow(theta, 0.8);

        double F_i = LineShapeFactor(f__ghz, OxygenData::f_0[i], delta_f__ghz, delta);

        // summation of terms...from Equation 2a, for oxygen
        N += S_i * F_i;
    }

    double N_D = NonresonantDebyeAttenuation(f__ghz, e__hPa, p__hPa, theta);

    double N_o = N + N_D;

    return N_o;
}

double OxygenRefractivitySimd(double f__ghz, double T__kelvin, double e__hPa, double p__hPa)
{
    double theta = 300 / T__kelvin;

    double N = 0;

    std::size_t i = 0;
    for (; i + doublev::size() <= OxygenData::f_0.size(); i += doublev::size())
    {
        // Equation 3, for oxygen
        const doublev a_1_simd = doublev(&OxygenData::a_1[i], stdx::element_aligned);
        const doublev a_2_simd = doublev(&OxygenData::a_2[i], stdx::element_aligned);
        const doublev S_i = a_1_simd * 1e-7 * p__hPa * pow(theta, 3) * exp(a_2_simd * (1 - theta));

        // compute the width of the line, Equation 6a, for oxygen
        const doublev a_3_simd = doublev(&OxygenData::a_3[i], stdx::element_aligned);
        const doublev a_4_simd = doublev(&OxygenData::a_4[i], stdx::element_aligned);
        doublev delta_f__ghz = a_3_simd * 1e-4 * (p__hPa * pow(theta, (0.8 - a_4_simd)) + 1.1 * e__hPa * theta);

        // modify the line width to account for Zeeman splitting of the oxygen lines
        // Equation 6b, for oxygen
        delta_f__ghz = sqrt(pow(delta_f__ghz, 2) + 2.25e-6);

        // correction factor due to interference effects in oxygen lines
        // Equation 7, for oxygen
        const doublev a_5_simd = doublev(&OxygenData::a_5[i], stdx::element_aligned);
        const doublev a_6_simd = doublev(&OxygenData::a_6[i], stdx::element_aligned);
        const doublev delta = (a_5_simd + a_6_simd * theta) * 1e-4 * (p__hPa + e__hPa) * pow(theta, 0.8);

        const doublev f_0_simd = doublev(&OxygenData::f_0[i], stdx::element_aligned);
        const doublev F_i = LineShapeFactorSimd(f__ghz, f_0_simd, delta_f__ghz, delta);

        // summation of terms...from Equation 2a, for oxygen
        N += stdx::reduce(S_i * F_i, std::plus{});
    }

    for (; i < OxygenData::f_0.size(); i++)
    {
        // Equation 3, for oxygen
        double S_i = OxygenData::a_1[i] * 1e-7 * p__hPa * pow(theta, 3) * exp(OxygenData::a_2[i] * (1 - theta));

        // compute the width of the line, Equation 6a, for oxygen
        double delta_f__ghz = OxygenData::a_3[i] * 1e-4 * (p__hPa * pow(theta, (0.8 - OxygenData::a_4[i])) + 1.1 * e__hPa * theta);

        // modify the line width to account for Zeeman splitting of the oxygen lines
        // Equation 6b, for oxygen
        delta_f__ghz = sqrt(pow(delta_f__ghz, 2) + 2.25e-6);

        // correction factor due to interference effects in oxygen lines
        // Equation 7, for oxygen
        double delta = (OxygenData::a_5[i] + OxygenData::a_6[i] * theta) * 1e-4 * (p__hPa + e__hPa) * pow(theta, 0.8);

        double F_i = LineShapeFactor(f__ghz, OxygenData::f_0[i], delta_f__ghz, delta);

        // summation of terms...from Equation 2a, for oxygen
        N += S_i * F_i;
    }

    double N_D = NonresonantDebyeAttenuation(f__ghz, e__hPa, p__hPa, theta);

    double N_o = N + N_D;

    return N_o;
}

/*=============================================================================
 |
 |  Description:  Imaginary part of the frequency-dependent complex
 |                refractivity due to water vapour.  See Equation (2b).
 |
 |        Input:  f__ghz        - Frequency, in GHz
 |                T__kelvin     - Temperature, in Kelvin
 |                e__hPa        - Water vapour partial pressure, in hPa
 |                p__hPa        - Dry air pressure, in hPa
 |
 |      Returns:  N_w           - Refractivity, in N-Units
 |
 *===========================================================================*/
double WaterVapourRefractivity(double f__ghz, double T__kelvin, double e__hPa, double P__hPa)
{
    double theta = 300 / T__kelvin;

    double N_w = 0;

    for (int i = 0; i < WaterVapourData::f_0.size(); i++)
    {
        // Equation 3, for water vapour
        double S_i = 0.1 * WaterVapourData::b_1[i] * e__hPa * pow(theta, 3.5) * exp(WaterVapourData::b_2[i] * (1 - theta));

        // compute the width of the line, Equation 6a, for water vapour
        double delta_f__ghz = 1e-4 * WaterVapourData::b_3[i] * (P__hPa * pow(theta, WaterVapourData::b_4[i]) + WaterVapourData::b_5[i] * e__hPa * pow(theta, WaterVapourData::b_6[i]));

        // modify the line width to account for Doppler broadening of water vapour lines
        // Equation 6b, for water vapour
        double term1 = 0.217 * pow(delta_f__ghz, 2) + (2.1316e-12 * pow(WaterVapourData::f_0[i], 2) / theta);
        delta_f__ghz = 0.535 * delta_f__ghz + sqrt(term1);

        // Equation 7, for water vapour
        double delta = 0;

        const double F_i = LineShapeFactor(f__ghz, WaterVapourData::f_0[i], delta_f__ghz, delta);

        // summation of terms...from Equation 2b, for water vapour
        N_w += S_i * F_i;
    }

    return N_w;
}

double WaterVapourRefractivitySimd(double f__ghz, double T__kelvin, double e__hPa, double P__hPa)
{
    double theta = 300 / T__kelvin;

    double N_w = 0;

    std::size_t i = 0;
    for (; i + doublev::size() <= WaterVapourData::f_0.size(); i += doublev::size())
    {
        // Equation 3, for water vapour
        const doublev b_1_simd = doublev(&WaterVapourData::b_1[i], stdx::element_aligned);
        const doublev b_2_simd = doublev(&WaterVapourData::b_2[i], stdx::element_aligned);
        const doublev S_i = 0.1 * b_1_simd * e__hPa * pow(theta, 3.5) * exp(b_2_simd * (1 - theta));

        // compute the width of the line, Equation 6a, for water vapour
        const doublev b_3_simd = doublev(&WaterVapourData::b_3[i], stdx::element_aligned);
        const doublev b_4_simd = doublev(&WaterVapourData::b_4[i], stdx::element_aligned);
        const doublev b_5_simd = doublev(&WaterVapourData::b_5[i], stdx::element_aligned);
        const doublev b_6_simd = doublev(&WaterVapourData::b_6[i], stdx::element_aligned);
        doublev delta_f__ghz = 1e-4 * b_3_simd * (P__hPa * pow(theta, b_4_simd) + b_5_simd * e__hPa * pow(theta, b_6_simd));

        // modify the line width to account for Doppler broadening of water vapour lines
        // Equation 6b, for water vapour
        const doublev f_0_simd = doublev(&WaterVapourData::f_0[i], stdx::element_aligned);
        const doublev term1 = 0.217 * pow(delta_f__ghz, 2) + (2.1316e-12 * pow(f_0_simd, 2) / theta);
        delta_f__ghz = 0.535 * delta_f__ghz + sqrt(term1);

        // Equation 7, for water vapour
        const doublev delta{0};

        const doublev F_i = LineShapeFactorSimd(f__ghz, f_0_simd, delta_f__ghz, delta);

        // summation of terms...from Equation 2b, for water vapour
        N_w += stdx::reduce(S_i * F_i, std::plus{});
    }

    for (; i < WaterVapourData::f_0.size(); i++)
    {
        // Equation 3, for water vapour
        double S_i = 0.1 * WaterVapourData::b_1[i] * e__hPa * pow(theta, 3.5) * exp(WaterVapourData::b_2[i] * (1 - theta));

        // compute the width of the line, Equation 6a, for water vapour
        double delta_f__ghz = 1e-4 * WaterVapourData::b_3[i] * (P__hPa * pow(theta, WaterVapourData::b_4[i]) + WaterVapourData::b_5[i] * e__hPa * pow(theta, WaterVapourData::b_6[i]));

        // modify the line width to account for Doppler broadening of water vapour lines
        // Equation 6b, for water vapour
        double term1 = 0.217 * pow(delta_f__ghz, 2) + (2.1316e-12 * pow(WaterVapourData::f_0[i], 2) / theta);
        delta_f__ghz = 0.535 * delta_f__ghz + sqrt(term1);

        // Equation 7, for water vapour
        double delta = 0;

        const double F_i = LineShapeFactor(f__ghz, WaterVapourData::f_0[i], delta_f__ghz, delta);

        // summation of terms...from Equation 2b, for water vapour
        N_w += S_i * F_i;
    }

    return N_w;
}
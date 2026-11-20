#include "p676.h"
#include <cmath>

namespace ITS {
namespace ITU {
namespace PSeries {
namespace p528 {

namespace {

/*=============================================================================
 |
 |  Description:  The exponents of theta that the lines of a gas use, without
 |                their repetitions. The line-by-line summation raises theta
 |                to an exponent of the line data for every line, but many
 |                lines share their exponents. Each power is computed once
 |                per summation and then looked up, which gives the very
 |                same values with a fraction of the calls of pow().
 |
 *===========================================================================*/
struct ExponentTable
{
    static const int MAX_EXPONENTS = 128;   // Size of the buffer for the powers

    vector<double> exponents;               // Distinct exponents
    vector<int> index_a;                    // For each line: its first exponent in exponents
    vector<int> index_b;                    // For each line: its second exponent in exponents
    bool usable = false;                    // false if the distinct exponents do not fit the buffer

    ExponentTable(const vector<double>& a, const vector<double>& b)
    {
        for (double value : a)
            index_a.push_back(Find(value));
        for (double value : b)
            index_b.push_back(Find(value));
        usable = (int)exponents.size() <= MAX_EXPONENTS;
    }

    int Find(double value)
    {
        for (size_t k = 0; k < exponents.size(); k++)
            if (exponents[k] == value)
                return (int)k;
        exponents.push_back(value);
        return (int)exponents.size() - 1;
    }

    // theta raised to each of the distinct exponents
    void Powers(double theta, double* powers) const
    {
        for (size_t k = 0; k < exponents.size(); k++)
            powers[k] = pow(theta, exponents[k]);
    }
};

// Exponents of the line width of oxygen, Equation 6a: 0.8 - a_4
const ExponentTable& OxygenExponents()
{
    static const ExponentTable table = [] {
        vector<double> width(OxygenData::a_4.size());
        for (size_t i = 0; i < width.size(); i++)
            width[i] = 0.8 - OxygenData::a_4[i];
        return ExponentTable(width, vector<double>());
    }();
    return table;
}

// Exponents of the line width of water vapour, Equation 6a: b_4 and b_6
const ExponentTable& WaterVapourExponents()
{
    static const ExponentTable table(WaterVapourData::b_4, WaterVapourData::b_6);
    return table;
}

// Squares of the line center frequencies of water vapour, for Equation 6b
const vector<double>& WaterVapourFrequenciesSquared()
{
    static const vector<double> squares = [] {
        vector<double> values(WaterVapourData::f_0.size());
        for (size_t i = 0; i < values.size(); i++)
            values[i] = pow(WaterVapourData::f_0[i], 2);
        return values;
    }();
    return squares;
}

} // namespace

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

    // powers of theta that do not depend on the line, or only through an
    // exponent that lines have in common
    const double theta_3 = pow(theta, 3);
    const double theta_08 = pow(theta, 0.8);
    const ExponentTable& exponents = OxygenExponents();
    double theta_width[ExponentTable::MAX_EXPONENTS];
    if (exponents.usable)
        exponents.Powers(theta, theta_width);

    for (int i = 0; i < OxygenData::f_0.size(); i++)
    {
        // Equation 3, for oxygen
        double S_i = OxygenData::a_1[i] * 1e-7 * p__hPa * theta_3 * exp(OxygenData::a_2[i] * (1 - theta));

        // compute the width of the line, Equation 6a, for oxygen
        double theta_a_4 = exponents.usable ? theta_width[exponents.index_a[i]] : pow(theta, (0.8 - OxygenData::a_4[i]));
        double delta_f__ghz = OxygenData::a_3[i] * 1e-4 * (p__hPa * theta_a_4 + 1.1  * e__hPa * theta);

        // modify the line width to account for Zeeman splitting of the oxygen lines
        // Equation 6b, for oxygen
        delta_f__ghz = sqrt(pow(delta_f__ghz, 2) + 2.25e-6);

        // correction factor due to interference effects in oxygen lines
        // Equation 7, for oxygen
        double delta = (OxygenData::a_5[i] + OxygenData::a_6[i] * theta) * 1e-4 * (p__hPa + e__hPa) * theta_08;

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

    // powers of theta that do not depend on the line, or only through an
    // exponent that lines have in common
    const double theta_35 = pow(theta, 3.5);
    const ExponentTable& exponents = WaterVapourExponents();
    double theta_width[ExponentTable::MAX_EXPONENTS];
    if (exponents.usable)
        exponents.Powers(theta, theta_width);
    const vector<double>& f_0_squared = WaterVapourFrequenciesSquared();

    for (int i = 0; i < WaterVapourData::f_0.size(); i++)
    {
        // Equation 3, for water vapour
        double S_i = 0.1 * WaterVapourData::b_1[i] * e__hPa * theta_35 * exp(WaterVapourData::b_2[i] * (1 - theta));

        // compute the width of the line, Equation 6a, for water vapour
        double theta_b_4 = exponents.usable ? theta_width[exponents.index_a[i]] : pow(theta, WaterVapourData::b_4[i]);
        double theta_b_6 = exponents.usable ? theta_width[exponents.index_b[i]] : pow(theta, WaterVapourData::b_6[i]);
        double delta_f__ghz = 1e-4 * WaterVapourData::b_3[i] * (P__hPa * theta_b_4 + WaterVapourData::b_5[i] * e__hPa * theta_b_6);

        // modify the line width to account for Doppler broadening of water vapour lines
        // Equation 6b, for water vapour
        double term1 = 0.217 * pow(delta_f__ghz, 2) + (2.1316e-12 * f_0_squared[i] / theta);
        delta_f__ghz = 0.535 * delta_f__ghz + sqrt(term1);

        // Equation 7, for water vapour
        double delta = 0;

        double F_i = LineShapeFactor(f__ghz, WaterVapourData::f_0[i], delta_f__ghz, delta);

        // summation of terms...from Equation 2b, for water vapour
        N_w += S_i * F_i;
    }

    return N_w;
}
}  // namespace p528
}  // namespace PSeries
}  // namespace ITU
}  // namespace ITS

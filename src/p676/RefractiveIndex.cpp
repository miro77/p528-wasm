#include "p676.h"
#include <cmath>

namespace ITS {
namespace ITU {
namespace PSeries {
namespace p528 {

/*=============================================================================
 |
 |  Description:  Compute the radio refractivity, n - 1. Unlike the
 |                refractive index it keeps its full precision, which
 |                matters when it is interpolated.
 |
 |        Input:  p__hPa        - Dry pressure, in hPa
 |                T__kelvin     - Temperature, in Kelvin
 |                e__hPa        - Water vapour pressure, in hPa
 |
 |      Returns:  n - 1         - Refractive index minus one
 |
 *===========================================================================*/
double RefractiveIndexMinusOne(double p__hPa, double T__kelvin, double e__hPa)
{
    // dry term of refractivity
    double N_dry = 77.6 * p__hPa / T__kelvin;

    // wet term of refractivity
    double N_wet = 72 * e__hPa / T__kelvin + 3.75e5 * e__hPa / pow(T__kelvin, 2);

    double N = N_dry + N_wet;

    return N * pow(10, -6);
}

/*=============================================================================
 |
 |  Description:  Compute the refractive index.
 |
 |        Input:  p__hPa        - Dry pressure, in hPa
 |                T__kelvin     - Temperature, in Kelvin
 |                e__hPa        - Water vapour pressure, in hPa
 |
 |      Returns:  n             - Refractive index
 |
 *===========================================================================*/
double RefractiveIndex(double p__hPa, double T__kelvin, double e__hPa)
{
    double n = 1 + RefractiveIndexMinusOne(p__hPa, T__kelvin, e__hPa);

    return n;
}

}  // namespace p528
}  // namespace PSeries
}  // namespace ITU
}  // namespace ITS

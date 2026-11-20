#include "p676.h"
#include "p835.h"
#include <cmath>

namespace ITS {
namespace ITU {
namespace PSeries {
namespace p528 {

namespace {

// The search for the grazing height stops when the two sides of Snell's law
// agree to this, see SlantPathAttenuation()
const double TERMINATE__KM = 0.001;

// A stored comparison could differ from the one in SlantPathAttenuation() in
// the last bit. Rays this close to the stopping criterion are left to
// SlantPathAttenuation(), so that both always stop at the same height.
const double MARGIN__KM = 1e-10;

// Deeper searches than this are left to SlantPathAttenuation() as well
const int MAX_DEPTH = 64;

} // namespace

/*=============================================================================
 |
 |  Description:  Cache for rays with a negative elevation angle between two
 |                fixed heights, see Section 2.2.2. Such a ray first
 |                descends from h_1 to a grazing height h_G and then rises
 |                to h_2.
 |
 |                SlantPathAttenuation() finds h_G by halving an interval
 |                until Snell's law holds to 1 m, and then traces from h_G
 |                to both heights. The result is not a smooth function of
 |                the elevation angle: h_G only takes the values of that
 |                search, and each of them moves the layers of both traces.
 |                So there is nothing to interpolate. But the heights the
 |                search visits form a binary tree that is the same for
 |                every elevation angle, and many rays of a profile end at
 |                the same h_G. The cache stores the tree: every height is
 |                evaluated once, and so are the two traces of every h_G.
 |                The results are those of SlantPathAttenuation().
 |
 |        Input:  f__ghz        - Frequency, in GHz
 |                h_1__km       - Height at which the rays start, in km
 |                h_2__km       - Height at which the rays end, in km
 |                table         - Layer properties for the frequency. It
 |                                has to outlive this object.
 |
 *===========================================================================*/
GrazingRayCache::GrazingRayCache(double f__ghz, double h_1__km, double h_2__km,
    SpecificAttenuationTable* table)
    : f__ghz_(f__ghz), h_1__km_(h_1__km), h_2__km_(h_2__km), table_(table)
{
    // the atmosphere of SlantPathAttenuation()
    config_.temperature = GlobalTemperature;
    config_.dry_pressure = GlobalPressure;
    config_.wet_pressure = GlobalWetPressure;

    // refractive index at h_1
    double p__hPa = config_.dry_pressure(h_1__km);
    double T__kelvin = config_.temperature(h_1__km);
    double e__hPa = config_.wet_pressure(h_1__km);

    n_1_ = RefractiveIndex(p__hPa, T__kelvin, e__hPa);
}

/*=============================================================================
 |
 |  Description:  Adds a height of the search to the tree.
 |
 |        Input:  h_G__km       - Candidate for the grazing height, in km
 |
 |      Returns:  Index of the new node
 |
 *===========================================================================*/
int GrazingRayCache::AddNode(double h_G__km)
{
    double p__hPa = config_.dry_pressure(h_G__km);
    double T__kelvin = config_.temperature(h_G__km);
    double e__hPa = config_.wet_pressure(h_G__km);

    double n_G = RefractiveIndex(p__hPa, T__kelvin, e__hPa);

    Node node;
    node.h_G__km = h_G__km;
    node.grazing_term = n_G * (a_0__km + h_G__km);
    nodes_.push_back(node);
    return (int)nodes_.size() - 1;
}

/*=============================================================================
 |
 |  Description:  Loss and length of the ray that starts at h_1 with a
 |                negative elevation angle and ends at h_2, like
 |                SlantPathAttenuation() with a table.
 |
 |        Input:  beta_1__rad   - Elevation angle (from zenith), in rad.
 |                                Has to be larger than 90 degrees.
 |
 |       Output:  A_gas__db     - Gaseous absorption along the ray, in dB
 |                a__km         - Ray length, in km
 |
 *===========================================================================*/
void GrazingRayCache::Lookup(double beta_1__rad, double* A_gas__db, double* a__km)
{
    const double start_term = n_1_ * (a_0__km + h_1__km_) * sin(beta_1__rad);

    // The search of SlantPathAttenuation(): it starts at the mid-point
    // between h_1 and the surface of the earth and halves its step.
    double delta = h_1__km_ / 2;
    if (nodes_.empty())
        AddNode(h_1__km_ - delta);

    int i_node = 0;
    for (int depth = 0; depth < MAX_DEPTH; depth++)
    {
        delta /= 2;

        double diff = nodes_[i_node].grazing_term - start_term;
        if (!(fabs(fabs(diff) - TERMINATE__KM) > MARGIN__KM))
            break;      // too close to call (or NaN)

        if (fabs(diff) <= TERMINATE__KM)
        {
            // converged on h_G. Trace in both directions with grazing angle
            if (!nodes_[i_node].traced)
            {
                double h_G__km = nodes_[i_node].h_G__km;
                SlantPathAttenuationResult result_1, result_2;
                double beta_graze__rad = PI / 2;
                RayTrace(f__ghz_, h_G__km, h_1__km_, beta_graze__rad, config_, &result_1, table_);
                RayTrace(f__ghz_, h_G__km, h_2__km_, beta_graze__rad, config_, &result_2, table_);

                nodes_[i_node].A_gas__db = result_1.A_gas__db + result_2.A_gas__db;
                nodes_[i_node].a__km = result_1.a__km + result_2.a__km;
                nodes_[i_node].traced = true;
            }
            *A_gas__db = nodes_[i_node].A_gas__db;
            *a__km = nodes_[i_node].a__km;
            return;
        }

        // next height of the search: below if the grazing term is too large
        int side = (diff > 0) ? 0 : 1;
        int i_child = nodes_[i_node].child[side];
        if (i_child < 0)
        {
            double h_G__km = nodes_[i_node].h_G__km;
            if (side == 0)
                h_G__km -= delta;
            else
                h_G__km += delta;
            i_child = AddNode(h_G__km);
            nodes_[i_node].child[side] = i_child;
        }
        i_node = i_child;
    }

    // not decided by the tree
    SlantPathAttenuationResult result;
    SlantPathAttenuation(f__ghz_, h_1__km_, h_2__km_, beta_1__rad, &result, table_);
    *A_gas__db = result.A_gas__db;
    *a__km = result.a__km;
}

}  // namespace p528
}  // namespace PSeries
}  // namespace ITU
}  // namespace ITS

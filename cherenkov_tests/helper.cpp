//
// Created by Matthew Dutson on 1/17/17.
//

#include "helper.h"
#include <TMath.h>

using namespace TMath;

namespace cherenkov_tests
{
    bool VectorsEqual(TVector3 actual, TVector3 expected, double fractional_err)
    {
        bool x_equal = ValuesEqual(actual.X(), expected.X(), fractional_err);
        bool y_equal = ValuesEqual(actual.Y(), expected.Y(), fractional_err);
        bool z_equal = ValuesEqual(actual.Z(), expected.Z(), fractional_err);
        return x_equal && y_equal && z_equal;
    }

    bool ValuesEqual(double actual, double expected, double fractional_err)
    {
        if (expected == 0)
        {
            return Abs(actual) < fractional_err;
        }
        else
        {
            return Abs((actual - expected) / expected) < fractional_err;
        }
    }
}
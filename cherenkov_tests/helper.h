//
// Created by Matthew Dutson on 1/17/17.
//

#ifndef test_helper_h
#define test_helper_h

#include <TVector3.h>

namespace cherenkov_tests
{
    /*
     * A function which will check whether two vectors are equal within acceptable error. The allowable fractional
     * between each component is specified.
     */
    bool VectorsEqual(TVector3 actual, TVector3 expected, double fractional_err);

    /*
     * Determines whether two decimals are equal within some acceptable fractional error.
     */
    bool ValuesEqual(double actual, double expected, double fractional_err);
}

#endif

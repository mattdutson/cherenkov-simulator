// Helper.h
//
// Author: Matthew Dutson
//
// Contains unit testing helper methods.

#ifndef HELPER_H
#define HELPER_H

#include <TVector3.h>

namespace cherenkov_simulator
{
    class Helper
    {
    public:

        /*
         * A function which will check whether two vectors are equal within acceptable error. The allowable fractional
         * between each component is specified.
         */
        static bool VectorsEqual(TVector3 actual, TVector3 expected, double fractional_err);

        /*
         * Determines whether two decimals are equal within some acceptable fractional error.
         */
        static bool ValuesEqual(double actual, double expected, double fractional_err);
    };
}

#endif

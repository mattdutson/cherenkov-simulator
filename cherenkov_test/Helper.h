// Helper.h
//
// Author: Matthew Dutson
//
// Definition of unit testing Helper class.

#ifndef HELPER_H
#define HELPER_H

#include <TVector3.h>

namespace cherenkov_simulator
{
    /*
     * Static class which contains unit testing helper methods.
     */
    class Helper
    {
    public:

        /*
         * A function which will check whether two vectors are equal within acceptable error. The allowable fractional
         * difference between each component is specified.
         */
        static bool VectorsEqual(TVector3 actual, TVector3 expected, double fractional_err);

        /*
         * Determines whether two decimals are equal within some acceptable fractional error.
         */
        static bool ValuesEqual(double actual, double expected, double fractional_err);
    };
}

#endif

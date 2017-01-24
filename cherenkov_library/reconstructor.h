// reconstructor.h
// cherenkov_library
//
// Created by Matthew Dutson on 1/8/17.
//

#ifndef reconstructor_h
#define reconstructor_h

#include "data_containers.h"
#include "geometric_objects.h"
#include "TRotation.h"
#include <istream>
#include "boost/property_tree/ptree.hpp"

namespace cherenkov_library
{
    class Reconstructor
    {
    public:

        void ParseFile(boost::property_tree::ptree config);

        TVector3 FitSDPlane(PhotonCount data);

    private:

        /*
         * A method to count the largest cluster of adjacent "true" values in a 2D array.
         */
        int LargestCluster(std::vector<std::vector<int>> not_counted);

        /*
         * A recursive method used by LargestCluster. Visits the item at the specified coordinates and any of its
         * neighbors. Returns the total number of visited elements.
         */
        int Visit(int i, int j, std::vector<std::vector<int>>* not_counted);

        // Parameters relating to the position and orientation of the detector relative to its surroundings
        Plane ground_plane;
        TRotation rotate_to_world;

        // Parameters defining the amount of night sky background noise (used to subtract the average noise level from
        // the signal).
        double sky_noise;
        double ground_noise;
    };
}

#endif

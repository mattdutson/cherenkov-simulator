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

        Reconstructor(boost::property_tree::ptree config);

        /*
         * Finds the shower-detector plane based on the distribution of data points.
         */
        Plane FitSDPlane(PhotonCount data);

        /*
         * Performs an ordinary monocular time profile reconstruction of the shower geometry. A ground impact point is
         * not used.
         */
        void TimeProfileFit(PhotonCount data, Plane sd_plane, double* t_0, double* impact_param, double* angle);

        /*
         * Apply triggering logic to the signal. Look for consecutive groups of pixels in each time bin which have
         * signals above some threshold. Also, elimiate any noise which is below some lower threshold. Return true if
         * any frames were triggered.
         */
        bool ApplyTriggering(PhotonCount* data);

    private:

        /*
         * A method to count the largest cluster of adjacent "true" values in a 2D array.
         */
        int LargestCluster(std::vector<std::vector<bool>> not_counted);

        /*
         * A recursive method used by LargestCluster. Visits the item at the specified coordinates and any of its
         * neighbors. Returns the total number of visited elements.
         */
        int Visit(int i, int j, std::vector<std::vector<bool>>* not_counted);

        // Parameters relating to the position and orientation of the detector relative to its surroundings
        Plane ground_plane;
        TRotation rotate_to_world;

        // Parameters defining the amount of night sky background noise
        double sky_noise;
        double ground_noise;

        // Parameters used when applying triggering logic
        double trigger_thresh;
        double hold_thresh;
        int trigger_clust;
    };
}

#endif

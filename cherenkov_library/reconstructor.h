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
#include <TGraph.h>
#include <TGraphErrors.h>

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
        TGraph MonocularFit(PhotonCount data, Plane sd_plane, double* t_0, double* impact_param, double* angle);

        /*
         * Performs a time profile reconstruction, but using the constraint of an impact point.
         */
        TGraph HybridFit(PhotonCount data, TVector3 impact_point, Plane sd_plane, double* t_0, double* impact_param,
                         double* angle);

        /*
         * Apply triggering logic to the signal. Look for consecutive groups of pixels in each time bin which have
         * signals above some threshold. Also, eliminate any noise which is below some lower threshold. Return true if
         * any frames were triggered.
         */
        bool ApplyTriggering(PhotonCount* data);

        /*
         * Performs the hybrid Cherenkov reconstruction. If the detector is not triggered, the "triggered" parameter is
         * set to false. If no valid ground impact point is found, then a standard reconstruction is done and
         * "ground_used" is set to false.
         */
        /*
         * Performs a standard time profile reconstruction. If the detector is not triggered, the "triggered" parameter
         * is set to false.
         */
        Shower Reconstruct(PhotonCount data, bool try_ground, bool* triggered, bool* ground_used);

    private:

        /*
         * Constructs the fit graph from data points.
         */
        TGraphErrors GetFitGraph(PhotonCount data, Plane sd_plane);

        /*
         * A method to count the largest cluster of adjacent "true" values in a 2D array.
         */
        int LargestCluster(std::vector<std::vector<bool>> not_counted);

        /*
         * A recursive method used by LargestCluster. Visits the item at the specified coordinates and any of its
         * neighbors. Returns the total number of visited elements.
         */
        int Visit(int i, int j, std::vector<std::vector<bool>>* not_counted);

        /*
         * Attempts to find the impact point of the shower. If this attempt fails, false is returned. Otherwise, true is
         * returned. We assume at this point that filters and triggering have been applied. Our condition is that some
         * pixel below the horizon must have seen a total number of photons which is more than three sigma from what we
         * would expect during that time frame.
         */
        bool FindGroundImpact(PhotonCount data, TVector3* impact);

        /*
         * Subtracts the average amount of noise from each pixel.
         */
        void SubtractNoise(PhotonCount* data);

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

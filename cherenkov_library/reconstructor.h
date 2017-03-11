// reconstructor.h
// cherenkov_library
//
// Created by Matthew Dutson on 1/8/17.
//

#ifndef reconstructor_h
#define reconstructor_h

#include <TRotation.h>
#include <istream>
#include <boost/property_tree/ptree.hpp>
#include <TGraph.h>
#include <TGraphErrors.h>

#include "data_containers.h"
#include "geometric_objects.h"

namespace cherenkov_library
{
    class Reconstructor
    {

    public:

        Reconstructor(boost::property_tree::ptree config);

        /*
         * Finds the shower-detector plane based on the distribution of data points. Returns a rotation to a frame in
         * which the shower-detector plane is the xy-plane, with the x-axis lying in the original xy-plane. This
         * rotation is assumed to start world frame, not the detector frame.
         */
        TRotation FitSDPlane(PhotonCount data);

        /*
         * Performs an ordinary monocular time profile reconstruction of the shower geometry. A ground impact point is
         * not used.
         */
        TGraphErrors
        MonocularFit(PhotonCount data, TRotation to_sd_plane, double* t_0, double* impact_param, double* angle);

        /*
         * Performs a time profile reconstruction, but using the constraint of an impact point.
         */
        TGraph
        HybridFit(PhotonCount data, TVector3 impact_point, TRotation to_sd_plane, double* t_0, double* impact_param,
                  double* angle);

        /*
         * Apply triggering logic to the signal. Look for consecutive groups of pixels in each time bin which have
         * signals above some threshold. Also, eliminate any noise which is below some lower threshold. Return true if
         * any frames were triggered.
         */
        bool ApplyTriggering(PhotonCount* data);

        /*
         * Determines whether the array of triggered tubes contains enough adjacent true values for the frame to be
         * triggered.
         */
        bool FrameTriggered(int t, std::vector<std::vector<std::vector<bool>>>* triggers);

        /*
         * Returns a 2D array of arrays, each of which contains true values for triggered tubes and false values for
         * untriggered tubes.
         */
        std::vector<std::vector<std::vector<bool>>> GetTriggeringMatrices(PhotonCount data);

        /*
         * Performs a reconstruction of the shower. If try_ground is true, then a hybrid Cherenkov reconstruction is
         * attempted. Otherwise, an ordinary monocular reconstruction is attempted. If the PhotonCount data wouldn't
         * result in the detector being triggered, then *triggered = false, and an empty shower is returned. If a valid
         * ground Cherenkov impact point is not found, then *ground_used = false, and an ordinary monocular
         * reconstruction is attempted.
         */
        Shower Reconstruct(PhotonCount data, bool try_ground, bool* triggered, bool* ground_used);

        /*
         * Subtracts the average amount of noise from each pixel.
         */
        void SubtractAverageNoise(PhotonCount* data);

        /*
         * Filters out any signal below three sigma. Assumes that the mean of the signal is zero (the average noise rate
         * has already been subtracted).
         */
        void ThreeSigmaFilter(PhotonCount* data);

        /*
         * Attempts to find the impact point of the shower. If this attempt fails, false is returned. Otherwise, true is
         * returned. We assume at this point that filters and triggering have been applied. Our condition is that some
         * pixel below the horizon must have seen a total number of photons which is more than three sigma from what we
         * would expect during that time frame.
         */
        bool FindGroundImpact(PhotonCount data, TVector3* impact);

    private:

        /*
         * Constructs the fit graph from data points.
         */
        TGraphErrors GetFitGraph(PhotonCount data, TRotation to_sd_plane);

        /*
         * A method to count the largest cluster of adjacent "true" values in a 2D array at a specified time bin.
         */
        int LargestCluster(int t, std::vector<std::vector<std::vector<bool>>>* not_counted);

        /*
         * A recursive method used by LargestCluster. Visits the item at the specified coordinates and any of its
         * neighbors. Returns the total number of visited elements.
         */
        int Visit(int i, int j, int t, std::vector<std::vector<std::vector<bool>>>* not_counted);

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

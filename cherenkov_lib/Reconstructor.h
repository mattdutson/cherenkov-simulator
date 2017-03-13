// Reconstructor.h
//
// Author: Matthew Dutson
//
// Contains a class used to reconstruct showers from photon arrival time data.

#ifndef RECONSTRUCTOR_H
#define RECONSTRUCTOR_H

#include <boost/property_tree/ptree.hpp>
#include <TRotation.h>
#include <TGraphErrors.h>

#include "DataStructures.h"
#include "Geometric.h"

namespace cherenkov_simulator
{
    class Reconstructor
    {
    public:

        /*
         * Constructs the Reconstructor from values in the configuration tree.
         */
        Reconstructor(boost::property_tree::ptree config);

        /*
         * Finds the shower-detector plane based on the distribution of data points. Returns a rotation to a frame in
         * which the shower-detector plane is the xy-plane, with the x-axis lying in the original xy-plane. This
         * rotation is assumed to start world frame, not the detector frame.
         */
        TRotation FitSDPlane(PhotonCount data);

        /*
         * Performs a reconstruction of the shower. If try_ground is true, then a hybrid Cherenkov reconstruction is
         * attempted. Otherwise, an ordinary monocular reconstruction is attempted. If the PhotonCount data wouldn't
         * result in the detector being triggered, then *triggered = false, and an empty shower is returned. If a valid
         * ground Cherenkov impact point is not found, then *ground_used = false, and an ordinary monocular
         * reconstruction is attempted.
         */
        Shower Reconstruct(PhotonCount data, bool try_ground, bool* triggered, bool* ground_used);

        /*
         * Performs an ordinary monocular time profile reconstruction of the shower geometry. A ground impact point is
         * not used.
         */
        TGraphErrors MonocularFit(PhotonCount data, TRotation to_sdp, double* t_0, double* r_p, double* psi);

        /*
         * Performs a time profile reconstruction, but using the constraint of an impact point.
         */
        TGraphErrors
        HybridFit(PhotonCount data, TVector3 impact, TRotation to_sdp, double* t_0, double* r_p, double* psi);

        /*
         * Attempts to find the impact point of the shower. If this attempt fails, false is returned. Otherwise, true is
         * returned. We assume at this point that filters and triggering have been applied. Our condition is that some
         * pixel below the horizon must have seen a total number of photons which is more than three sigma from what we
         * would expect during that time frame.
         */
        bool FindGroundImpact(PhotonCount data, TVector3* impact);

        /*
         * Adds Poisson-distributed background noise to the signal.
         */
        void AddNoise(PhotonCount* data);

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
         * Apply triggering logic to the signal. Look for consecutive groups of pixels in each time bin which have
         * signals above some threshold. Also, eliminate any noise which is below some lower threshold. Return true if
         * any frames were triggered.
         */
        bool ApplyTriggering(PhotonCount* data);

    private:

        friend class ReconstructorTest;

        // Levels of night sky background noise - cgs, sr
        constexpr static double global_sky_noise = 4.924e6;
        constexpr static double global_ground_noise = 4.924e5;

        // Parameters relating to the position and orientation of the detector relative to its surroundings - cgs
        Plane ground_plane;
        TRotation rotate_to_world;

        // Detector-specific levels of night sky background noise - cgs, sr
        double sky_noise;
        double ground_noise;

        // Parameters used when applying triggering logic
        double trigger_thresh;
        double hold_thresh;
        int trigger_clust;

        // A general-purpose random number generator
        TRandom3 rng;

        /*
         * Constructs the fit graph from data points.
         */
        TGraphErrors GetFitGraph(PhotonCount data, TRotation to_sdp);

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
         * A method to count the largest cluster of adjacent "true" values in a 2D array at a specified time bin.
         */
        int LargestCluster(int t, std::vector<std::vector<std::vector<bool>>>* not_counted);

        /*
         * A recursive method used by LargestCluster. Visits the item at the specified coordinates and any of its
         * neighbors. Returns the total number of visited elements.
         */
        int Visit(int i, int j, int t, std::vector<std::vector<std::vector<bool>>>* not_counted);
    };
}

#endif

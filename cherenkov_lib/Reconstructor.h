// Reconstructor.h
//
// Author: Matthew Dutson
//
// Contains a class used to reconstruct showers from photon arrival time data.

#ifndef RECONSTRUCTOR_H
#define RECONSTRUCTOR_H

#include <queue>
#include <boost/property_tree/ptree.hpp>
#include <TRotation.h>
#include <TGraphErrors.h>
#include <TMatrixDSym.h>

#include "Utility.h"
#include "DataStructures.h"
#include "Geometric.h"

namespace cherenkov_simulator
{
    class Reconstructor
    {
    public:

        struct Result
        {
            bool trigger;
            bool impact;
            Shower mono;
            Shower ckv;

            /*
             * Creates a header for rows of data created with ToString().
             */
            static std::string Header();

            /*
             * Creates a string with comma separated fields, with Shower represented by Shower.ToString().
             */
            std::string ToString();
        };

        /*
         * Constructs the Reconstructor from values in the configuration tree.
         */
        Reconstructor(const boost::property_tree::ptree& config);

        /*
         * Performs both a monocular and Cherenkov reconstruction, storing output in a Result data structure. If the
         * detector was not triggered, Result.triggered = false. If there was not visible impact point,
         * Result.cherenkov = false.
         */
        Result Reconstruct(PhotonCount& data);

        /*
         * Adds Poisson-distributed background noise to the signal.
         */
        void AddNoise(PhotonCount& data);

        /*
         * Attempts to isolate signal from noise by subtracting the background level, applying triggering, removing
         * anything below three sigma, and performing a recursive search from triggered pixels.
         */
        void ClearNoise(PhotonCount& data);

    private:

        friend class ReconstructorTest;

        // Levels of night sky background noise - cgs, sr
        constexpr static double global_sky_noise = 4.924e6;
        constexpr static double global_ground_noise = 4.924e5;

        // Parameters relating to the position and orientation of the detector relative to its surroundings - cgs
        Plane ground;
        TRotation to_world;

        // Detector-specific levels of night sky background noise - cgs, sr
        double sky_noise;
        double ground_noise;

        // Parameters used when applying triggering logic and noise reduction
        double trigger_thresh;
        double noise_thresh;
        int trigger_clust;
        double impact_buffer;
        double plane_dev;

        // A general-purpose random number generator
        TRandom3 rng;

        /*
         * Performs an ordinary monocular time profile reconstruction of the shower geometry. A ground impact point is
         * not used.
         */
        Shower MonocularFit(const PhotonCount& data, TRotation to_sdp, std::string graph_file = "");

        /*
         * Performs a time profile reconstruction, but using the constraint of an impact point.
         */
        Shower HybridFit(const PhotonCount& data, TVector3 impact, TRotation to_sdp, std::string graph_file = "");

        /*
         * Finds the shower-detector plane based on the distribution of data points. Returns a rotation to a frame in
         * which the shower-detector plane is the xy-plane, with the x-axis lying in the original xy-plane. This
         * rotation is assumed to start world frame, not the detector frame.
         */
        TRotation FitSDPlane(const PhotonCount& data, const Bool3D* mask = nullptr);

        /*
         * Finds the eigenvector of the symmetric matrix with the smallest eigenvalue.
         */
        TVector3 MinValVec(TMatrixDSym matrix);

        /*
         * Attempts to find the impact point of the shower. If this attempt fails, false is returned. Otherwise, true is
         * returned. We assume at this point that filters and triggering have been applied. Our condition is that some
         * pixel below the horizon must have seen a total number of photons which is more than three sigma from what we
         * would expect during that time frame.
         */
        bool FindGroundImpact(PhotonCount& data, TVector3& impact);

        /*
         * Constructs the fit graph from data points.
         */
        TGraphErrors GetFitGraph(const PhotonCount& data, TRotation to_sdp);

        /*
         * Subtracts the average amount of noise from each pixel.
         */
        void SubtractAverageNoise(PhotonCount& data);

        /*
         * Apply triggering logic to the signal. Look for consecutive groups of pixels in each time bin which have
         * signals above some threshold. Also, eliminate any noise which is below some lower threshold. Return true if
         * any frames were triggered.
         */
        std::vector<bool> GetTriggeringState(PhotonCount& data);

        /*
         * Modify the set of triggered pixels/times to contain the subset of triggered pixels/times which are within
         * some angle of an estimated shower-detector plane.
         */
        void FindPlaneSubset(const PhotonCount& data, Bool3D& triggered);

        /*
         * Determines whether the input direction is near enough to the plane. The maximum angular deviation from the
         * plane is defined in the config.
         */
        bool NearPlane(TRotation to_plane, TVector3 direction);

        /*
         * Determines whether the detector was triggered by iterating through trig_state and determining if there are
         * any "true" values.
         */
        bool DetectorTriggered(PhotonCount& data);

        /*
         * Returns a 2D array of arrays, each of which contains true values for tubes above the specified multiple of
         * sigma, and false values for all those below.
         */
        Bool3D GetThresholdMatrices(PhotonCount& data, double sigma_mult, bool use_below_horiz = true);

        /*
         * Constructs a shower based on the results of the time profile reconstruction.
         */
        Shower MakeShower(double t_0, double r_p, double psi, TRotation to_sd_plane);

        /*
         * Checks whether the specified point lies within the field of minus some buffer defined in the config. Assumes
         * the direction is in the detector frame.
         */
        bool PointWithinView(TVector3 direction, const PhotonCount& data);

        void VisitSpaceAdj(ULong x, ULong y, ULong t, std::queue<std::array<ULong, 3>>& front, const Bool3D& good);

        void VisitTimeAdj(ULong x, ULong y, ULong t, std::queue<std::array<ULong, 3>>& front, const Bool3D& good);

        void SafePush(ULong x, ULong y, ULong t, std::queue<std::array<ULong, 3>>& front, const Bool3D& good);
    };
}

#endif

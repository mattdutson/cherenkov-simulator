// Reconstructor.h
//
// Author: Matthew Dutson
//
// Contains a class used to reconstruct showers from photon arrival time data.

#ifndef RECONSTRUCTOR_H
#define RECONSTRUCTOR_H

#include <queue>
#include <array>
#include <boost/property_tree/ptree.hpp>
#include <TRotation.h>
#include <TGraphErrors.h>
#include <TMatrixDSym.h>

#include "DataStructures.h"
#include "Geometric.h"
#include "Utility.h"

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
            std::string ToString() const;
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
        Result Reconstruct(const PhotonCount& data, const Bool1D& trig) const;

        /*
         * Adds Poisson-distributed background noise to the signal.
         */
        void AddNoise(PhotonCount& data) const;

        /*
         * Attempts to isolate signal from noise by subtracting the background level, applying triggering, removing
         * anything below three sigma, and performing a recursive search from triggered pixels.
         */
        Bool1D ClearNoise(PhotonCount& data) const;

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

        /*
         * Performs an ordinary monocular time profile reconstruction of the shower geometry. A ground impact point is
         * not used.
         */
        Shower MonocularFit(const PhotonCount& data, TRotation to_sdp, std::string graph_file = "") const;

        /*
         * Performs a time profile reconstruction, but using the constraint of an impact point.
         */
        Shower HybridFit(const PhotonCount& data, TVector3 impact, TRotation to_sdp, std::string graph_file = "") const;

        /*
         * Finds the shower-detector plane based on the distribution of data points. Returns a rotation to a frame in
         * which the shower-detector plane is the xy-plane, with the x-axis lying in the original xy-plane. This
         * rotation is assumed to start world frame, not the detector frame.
         */
        TRotation FitSDPlane(const PhotonCount& data, const Bool3D* mask = nullptr) const;

        /*
         * Finds the eigenvector of the symmetric matrix with the smallest eigenvalue.
         */
        TVector3 MinValVec(TMatrixDSym matrix) const;

        /*
         * Attempts to find the impact point of the shower. If this attempt fails, false is returned. Otherwise, true is
         * returned. We assume at this point that filters and triggering have been applied. Our condition is that some
         * pixel below the horizon must have seen a total number of photons which is more than three sigma from what we
         * would expect during that time frame.
         */
        bool FindGroundImpact(const PhotonCount& data, TVector3& impact) const;

        /*
         * Constructs the fit graph from data points.
         */
        TGraphErrors GetFitGraph(const PhotonCount& data, TRotation to_sdp) const;

        /*
         * Subtracts the average amount of noise from each pixel.
         */
        void SubtractAverageNoise(PhotonCount& data) const;

        /*
         * Apply triggering logic to the signal. Look for consecutive groups of pixels in each time bin which have
         * signals above some threshold. Also, eliminate any noise which is below some lower threshold. Return true if
         * any frames were triggered.
         */
        Bool1D GetTriggeringState(const PhotonCount& data) const;

        /*
         * Visits all spatially adjacent pixels to the (x, y, t) point passed, pushing them to the queue. They are also
         * marked as visited in the not_visited structure.
         */
        static void VisitSpaceAdj(size_t x, size_t y, size_t t, std::list<std::array<size_t, 3>>& front, Bool3D& not_visited);

        /*
         * Visits all spatially adjacent temporally to the (x, y, t) point passed, pushing them to the queue. They are
         * also marked as visited in the not_visited structure.
        */
        static void VisitTimeAdj(size_t x, size_t y, size_t t, std::list<std::array<size_t, 3>>& front, Bool3D& not_visited);

        /*
         * Pushes the specified (x, y, z) point to the queue, first checking that the point lies within appropriate
         * bounds and has not yet been visited. not_visited[x][y][t] is set to false.
         */
        static void VisitPush(size_t x, size_t y, size_t t, std::list<std::array<size_t, 3>>& front, Bool3D& not_visited);

        /*
         * Modify the set of triggered pixels/times to contain the subset of triggered pixels/times which are within
         * some angle of an estimated shower-detector plane.
         */
        void FindPlaneSubset(const PhotonCount& data, Bool3D& triggered) const;

        /*
         * Determines whether the input direction is near enough to the plane. The maximum angular deviation from the
         * plane is defined in the config.
         */
        bool NearPlane(TRotation to_plane, TVector3 direction) const;

        /*
         * Determines whether the detector was triggered by iterating through trig_state and determining if there are
         * any "true" values.
         */
        bool DetectorTriggered(const Bool1D& trig_state) const;

        /*
         * Returns a 2D array of arrays, each of which contains true values for tubes above the specified multiple of
         * sigma, and false values for all those below.
         */
        Bool3D GetThresholdMatrices(const PhotonCount& data, double sigma_mult, bool use_below_horiz = true) const;

        /*
         * Constructs a shower based on the results of the time profile reconstruction.
         */
        static Shower MakeShower(double t_0, double r_p, double psi, TRotation to_sd_plane);
    };
}

#endif

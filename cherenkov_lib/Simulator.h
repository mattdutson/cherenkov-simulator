// Simulator.h
//
// Author: Matthew Dutson
//
// Contains the definition of the Simulator class, which performs shower simulation and ray tracing.

#ifndef simulator_h
#define simulator_h

#include <boost/property_tree/ptree.hpp>
#include <TF1.h>
#include <TRandom3.h>
#include <TRotation.h>

#include "DataStructures.h"
#include "Geometric.h"

namespace cherenkov_lib
{
    /*
     * A class which performs the majority of the shower simulation and contains most simulation parameters.
     */
    class Simulator
    {
    public:

        /*
         * Takes a parsed XML object and attempts to extract required simulation parameters. If this fails, an
         * exception is thrown which specifies the name of the missing parameter.
         */
        Simulator(boost::property_tree::ptree config);

        /*
         * Simulate the motion of the shower from its current point to the ground, emitting fluorescence and Cherenkov
         * photons at each depth step. Ray trace these photons through the Schmidt detector and record their impact
         * positions. Add Poisson-distributed background noise.
         */
        PhotonCount SimulateShower(Shower shower);

        /*
         * Adds Poisson-distributed background noise to the signal.
         */
        void AddNoise(PhotonCount* photon_count);

    private:

        class CherenkovFunc
        {
        public:

            CherenkovFunc(Simulator* sim);

            double operator()(double* x, double* p);

        private:

            Simulator* sim;

        };

        friend class CherenkovFunc;

        /*
         * Simulate the production and detection of the fluorescence photons.
         */
        void ViewFluorescencePhotons(Shower shower, PhotonCount* photon_count);

        /*
         * Simulate the production and detection of the Cherenkov photons. Only Cherenkov photons reflected from the
         * ground are recorded (no back scattering).
         */
        void ViewCherenkovPhotons(Shower shower, Plane ground_plane, PhotonCount* photon_count);

        /*
         * Takes a photon which is assumed to lie at the corrector plate and simulates its motion through the detector
         * optics. If the photon is somehow blocked or doesn't reach the photomultiplier array, no change to the photon
         * count structure is made. Otherwise, the appropriate bin of the photon counter is incremented. Takes a
         * parameter which represents the rate of computational thinning. This is passed to the photon count container
         * to allow it to increment bins by the correct amount.
         */
        void SimulateOptics(Ray photon, PhotonCount* photon_count, double thinning);

        /*
         * Generates a random point on the circle of the refracting lens.
         */
        TVector3 RandomStopImpact();

        /*
         * Refracts a ray across the Schmidt corrector. The Schmidt corrector is assumed to have zero thickness.
         */
        bool DeflectFromLens(Ray* photon);

        /*
         * Takes a ray which has just been refracted by the corrector. Finds the point on the mirror where that ray will
         * reflect. If the ray misses the mirror, false is returned.
         */
        bool MirrorImpactPoint(Ray ray, TVector3* point);

        /*
         * Returns the normal vector at some point on the mirror. Behavior is undefined if the passed point is not on or
         * near the mirror.
         */
        TVector3 MirrorNormal(TVector3 point);

        /*
         * Finds the point where some ray will impact the photomultiplier surface. This can be used both to check
         * whether photons are blocked by the photomultipliers and to find where they are detected after being reflected
         * by the mirror. Returns false if the ray will not hit the photomultiplier array.
         */
        bool CameraImpactPoint(Ray ray, TVector3* point);

        /*
         * Finds the points where a ray will or has intersected with a sphere centered at the origin. If the ray does
         * not intersect with the sphere, "point" is set to (0, 0, 0) and false is returned. Otherwise, "point" is set
         * to the intersection with the smallest (negative) z-coordinate and "true" is returned.
         */
        bool BackOriginSphereImpact(Ray ray, TVector3* point, double radius);

        /*
         * Determines the total number of Fluorescence photons produced by the shower at a particular point. Takes the
         * distance traveled by the shower due to the fact that the form for fluorescence yield (Stratton 4.2) gives the
         * number of photons per electron per unit length.
         *
         * EDIT: No longer takes the distance due to the fact that a modified version of the fluorescence yield formula
         * is being used. See notes on unit consolidation.
         */
        int NumberFluorescencePhotons(Shower shower);

        /*
         * Determines the total number of Cherenkov photons produced by the shower at a particular point. This doesn't
         * need the distance traveled because the form for Cherenkov yield gives the number of photons per electron per
         * slant depth.
         */
        int NumberCherenkovPhotons(Shower shower);

        /*
         * Calculates the number of charged particles in the shower using the Gaiser-Hilles profile.
         */
        double GaiserHilles(Shower shower);

        /*
         * Calculates the effective ionization loss rate for a shower (alpha_eff).
         */
        double IonizationLossRate(Shower shower);

        /*
         * Calculates the Cherenkov threshold energy at some height. Uses the AtmosphereDelta function.
         */
        double EThresh(Shower shower);

        /*
         * Calculates how large, as a fraction of a sphere, the detector stop appears from some point. This accounts
         * both for the inverse square dependance and the orientation of the detector.
         */
        double SphereFraction(TVector3 view_point);

        /*
         * Returns the product of the quantum efficiency, filter transmittance, and mirror reflectance.
         */
        double DetectorEfficiency();

        /*
         * Creates a Cherenkov photon with a randomly-assigned direction (the direction follows a e^-theta/sin(theta)
         * distribution.
         */
        Ray GenerateCherenkovPhoton(Shower shower);

        /*
         * Calculates the critical angle in the expression for the Cherenkov angular distribution.
         */
        double ThetaC(Shower shower);

        /*
         * Generates a time which is randomly offset from the shower time.
         */
        double JitteredTime(Shower shower);

        // Parameters related to the behavior of the simulation
        double depth_step;
        double time_bin;
        double fluor_thin;
        double ckv_thin;

        // Parameters relating to the position and orientation of the detector relative to its surroundings
        Plane ground_plane;
        TRotation rotate_to_world;

        // The constant temperature of the atmosphere
        double atmosphere_temp;

        // Parameters defining properties of the detector optics
        double refrac_lens;
        double mirror_radius;
        double stop_diameter;
        double mirror_size;
        double cluster_diameter;
        int n_pmt_across;

        // Parameters in the GH profile
        double gh_lambda;

        // Parameters used when calculating the fluorescence yield
        double fluor_a1;
        double fluor_a2;
        double fluor_b1;
        double fluor_b2;
        double dep_1_4;

        // Parameters used when calculating the effective ionization loss rate
        double ion_c1;
        double ion_c2;
        double ion_c3;
        double ion_c4;
        double ion_c5;

        // Parameters used when calculating theta_c in the Cherenkov angular distribution
        double ckv_k1;
        double ckv_k2;

        // Parameters used when calculating the Cherenkov yield
        double lambda_min;
        double lambda_max;

        // Parameters in the electron energy spectrum
        double fe_a11;
        double fe_a12;
        double fe_a21;
        double fe_a22;
        double fe_k0;
        double fe_k1;
        double fe_k2;

        // Physics constants
        double mass_e;
        double fine_struct;

        // Parameters defining the amount of night sky background noise
        double sky_noise;
        double ground_noise;

        // Parameters which describe inefficiencies in the equipment
        double mirror_reflect;
        double filter_transmit;
        double quantum_eff;

        // A general purpose random number generator
        TRandom3 rng;

        // The Cherenkov functor and corresponding TF1 for integration
        CherenkovFunc ckv_func;
        TF1 ckv_integrator;
    };
}

#endif

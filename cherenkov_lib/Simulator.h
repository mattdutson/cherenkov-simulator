// Simulator.h
//
// Author: Matthew Dutson
//
// Contains the definition of the Simulator class, which performs shower simulation and ray tracing.

#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <boost/property_tree/ptree.hpp>
#include <TF1.h>
#include <TRandom3.h>
#include <TRotation.h>

#include "DataStructures.h"
#include "Geometric.h"
#include "Utility.h"

namespace cherenkov_simulator
{
    /*
     * A class which performs the majority of the shower simulation and contains most simulation parameters.
     */
    class Simulator
    {
    public:

        /*
         * Constructs the Simulator from values in the configuration tree.
         */
        Simulator(const boost::property_tree::ptree& config);

        /*
         * Simulate the motion of the shower from its current point to the ground, emitting fluorescence and Cherenkov
         * photons at each depth step. Ray trace these photons through the Schmidt detector and record their impact
         * positions. Add Poisson-distributed background noise.
         */
        PhotonCount SimulateShower(Shower shower) const;

    private:

        friend class SimulatorTest;

        /*
         * Represents the integrand of the Cherenkov yield.
         */
        class CherenkovFunc
        {
        public:

            /*
             * Defined as required by TF1.Integral(). Returns the value of the integrand.
             * p[0] = age, p[1] = rho, p[2] = delta
             */
            double operator()(double* x, double* p);

        private:

            // Parameters in the Cherenkov yield - cgs
            constexpr static double lambda_min = 3.0e-5;
            constexpr static double lambda_max = 4.0e-5;

            // Parameters in the electron energy spectrum - MeV
            constexpr static double fe_a11 = 6.42522;
            constexpr static double fe_a12 = 1.53183;
            constexpr static double fe_a21 = 168.168;
            constexpr static double fe_a22 = 42.1368;
            constexpr static double fe_k0 = 1.45098e-1;
            constexpr static double fe_k1 = 6.20114;
            constexpr static double fe_k2 = -5.96851e-1;
        };

        // Miscellaneous parameters - K
        constexpr static double atm_temp = 273.0;
        constexpr static double refrac_lens = 1.52;

        // Parameters in the fluorescence yield - MeV, cgs, K
        constexpr static double fluor_a1 = 890.0;
        constexpr static double fluor_a2 = 550.0;
        constexpr static double fluor_b1 = 1850.0;
        constexpr static double fluor_b2 = 6500.0;
        constexpr static double dep_1_4 = 1.6;

        // Parameters in the effective ionization loss rate - MeV, cgs
        constexpr static double ion_c1 = 3.90883;
        constexpr static double ion_c2 = 1.05301;
        constexpr static double ion_c3 = 9.91717;
        constexpr static double ion_c4 = 2.41715;
        constexpr static double ion_c5 = 0.13180;

        // Parameters used when calculating theta_c in the Cherenkov angular distribution - MeV
        constexpr static double ckv_k1 = 0.83;
        constexpr static double ckv_k2 = -0.67;

        // Parameters which describe inefficiencies in the equipment
        constexpr static double mirror_reflect = 0.80; // 1.0; // 0.80
        constexpr static double filter_transmit = 1.0;
        constexpr static double quantum_eff = 0.15; // 1.0; // 0.15

        // Parameters related to the behavior of the simulation - cgs
        double depth_step;
        double fluor_thin;
        double ckv_thin;
        double back_toler;

        // Miscellaneous non-constant parameters
        Plane ground_plane;
        TRotation rotate_to_world;
        CherenkovFunc ckv_func;
        TF1 ckv_integrator;
        PhotonCount::Params count_params;

        // Setup of the detector - cgs
        double mirror_radius;
        double stop_diameter;
        double mirror_size;
        double cluster_diameter;

        /*
         * Simulate the production and detection of the fluorescence photons.
         */
        void ViewFluorescencePhotons(Shower shower, PhotonCount& photon_count) const;

        /*
         * Simulate the production and detection of the Cherenkov photons. Only Cherenkov photons reflected from the
         * ground are recorded (no back scattering).
         */
        void ViewCherenkovPhotons(Shower shower, Plane ground_plane, PhotonCount& photon_count, TF1 integrator) const;

        /*
         * Determines the total number of Fluorescence photons produced by the shower at a particular point. Takes the
         * distance traveled by the shower due to the fact that the form for fluorescence yield (Stratton 4.2) gives the
         * number of photons per electron per unit length.
         */
        int NumberFluorescencePhotons(Shower shower) const;

        /*
         * Determines the total number of Cherenkov photons produced by the shower at a particular point. This doesn't
         * need the distance traveled because the form for Cherenkov yield gives the number of photons per electron per
         * slant depth.
         */
        int NumberCherenkovPhotons(Shower shower, TF1 integrator) const;

        /*
         * Takes a photon which is assumed to lie at the corrector plate and simulates its motion through the detector
         * optics. If the photon is somehow blocked or doesn't reach the photomultiplier array, no change to the photon
         * count structure is made. Otherwise, the appropriate bin of the photon counter is incremented. Takes a
         * parameter which represents the rate of computational thinning. This is passed to the photon count container
         * to allow it to increment bins by the correct amount.
         */
        void SimulateOptics(Ray photon, PhotonCount& photon_count, double thinning) const;

        /*
         * Generates a random point on the circle of the refracting lens.
         */
        TVector3 RandomStopImpact() const;

        /*
         * Refracts a ray across the Schmidt corrector. The Schmidt corrector is assumed to have zero thickness.
         */
        bool DeflectFromLens(Ray& photon) const;

        /*
         * Takes a ray which has just been refracted by the corrector. Finds the point on the mirror where that ray will
         * reflect. If the ray misses the mirror, false is returned.
         */
        bool MirrorImpactPoint(Ray ray, TVector3& point) const;

        /*
         * Returns the normal vector at some point on the mirror. Behavior is undefined if the passed point is not on or
         * near the mirror.
         */
        TVector3 MirrorNormal(TVector3 point) const;

        /*
         * Finds the point where some ray will impact the photomultiplier surface. This can be used both to check
         * whether photons are blocked by the photomultipliers and to find where they are detected after being reflected
         * by the mirror. Returns false if the ray will not hit the photomultiplier array.
         */
        bool CameraImpactPoint(Ray ray, TVector3& point) const;

        /*
         * Calculates the effective ionization loss rate for a shower (alpha_eff).
         */
        double IonizationLossRate(Shower shower) const;

        /*
         * Calculates how large, as a fraction of a sphere, the detector stop appears from some point. This accounts
         * both for the inverse square dependance and the orientation of the detector.
         */
        double SphereFraction(TVector3 view_point) const;

        /*
         * Returns the product of the quantum efficiency, filter transmittance, and mirror reflectance.
         */
        double DetectorEfficiency() const;

        /*
         * Creates a Cherenkov photon with a randomly-assigned direction (the direction follows a e^-theta/sin(theta)
         * distribution.
         */
        Ray GenerateCherenkovPhoton(Shower shower) const;

        /*
         * Calculates the critical angle in the expression for the Cherenkov angular distribution.
         */
        double ThetaC(Shower shower) const;

        /*
         * Generates a time which is randomly offset from the shower time.
         */
        Ray JitteredRay(Shower shower, TVector3 direction) const;

        /*
         * Determines the time when we want to start recording photons for the shower. This is calculated by taking the
         * time it would take for a photon to travel directly from the shower's current position to the detector.
         */
        double MinTime(Shower shower) const;

        /*
         * Determines the time when we want to stop recording photons from the shower. This is calculated by the time it
         * would take for a photon to travel along the shower axis to the ground plus some multiple of the time it would
         * take for that photon to then reach the detector.
         */
        double MaxTime(Shower shower) const;

        /*
         * Finds the points where a ray will or has intersected with a sphere centered at the origin. If the ray does
         * not intersect with the sphere, "point" is set to (0, 0, 0) and false is returned. Otherwise, "point" is set
         * to the intersection with the smallest (negative) z-coordinate and "true" is returned.
         */
        static bool NegSphereImpact(Ray ray, TVector3& point, double radius);
    };
}

#endif

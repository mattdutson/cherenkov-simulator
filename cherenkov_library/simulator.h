// simulator.h
// cherenkov_simulator
//
// Created by Matthew Dutson on 9/8/16.
//
//

#ifndef simulator_h
#define simulator_h

#include "data_containers.h"
#include "geometric_objects.h"
#include "utility.h"
#include "TF1.h"
#include "Math/Transform3D.h"
#include "TRandom3.h"
#include "TRotation.h"
#include <istream>

namespace cherenkov_simulator
{
    class Simulator
    {

    private:

        // Parameters related to the behavior of the simulation
        double depth_step;
        double time_bin;

        // Parameters relating to the position and orientation of the detector relative to its surroundings
        Plane ground_plane;
        TRotation rotate_to_world;

        // Parameters defining properties of the atmosphere
        double scale_height;
        double rho_0;
        double delta_0;

        // Parameters defining properties of the detector optics
        double refrac_lens;
        double mirror_radius;
        double stop_size;
        double mirror_size;
        double cluster_size;
        int n_pmt_across;

        // Parameters in the GH profile
        double gh_lambda;

        // Parameters used when calculating the fluorescence yield
        double fluor_a1, fluor_a2;
        double fluor_b1, fluor_b2;

        // Parameters used when calculating the effective ionization loss rate
        double ion_c1, ion_c2, ion_c3, ion_c4, ion_c5;

        // Parameters used when calculating theta_c in the Cherenkov angular distribution
        double ckv_k1, ckv_k2;

        // Parameters used when calculating the Cherenkov yield
        double lambda_min, lambda_max;

        // Parameters in the electron energy spectrum
        double fe_a11, fe_a12;
        double fe_a21, fe_a22;
        double fe_k0, fe_k1, fe_k2;

        // Physics constants
        double mass_e;
        double fine_struct;

        // A general purpose random number generator
        TRandom3 rng;

        void ViewFluorescencePhotons(Shower shower, PhotonCount* photon_count);

        void ViewCherenkovPhotons(Shower shower, Plane ground_plane, PhotonCount* photon_count);

        void SimulateOptics(Ray photon, PhotonCount* photon_count);

        /*
         * Finds the amount that the position of a given shower should be incremented based on the local atmospheric
         * density.
         */
        double FindDistance(Shower shower);

        /*
         * Adds Poisson-distributed background noise to the signal.
         */
        void AddNoise(PhotonCount* photon_count);

        TVector3 RandomStopImpact();

        TVector3 MirrorNormal(TVector3 point);

        bool MirrorImpactPoint(Ray ray, TVector3* point);

        bool CameraImpactPoint(Ray ray, TVector3* point);

        void DeflectFromLens(Ray* photon);

        /*
         * Calculates the intensity of the shower using the Gaiser-Hilles profile.
         */
        double GaiserHilles(Shower shower);

        /*
         * Calculates the fraction of photons radiated at a certain point which would be captured by our detector. The
         * point must be in the world reference frame (detector at origin with +z perpendicular to earth surface).
         */
        double PhotonFraction(TVector3 view_point);

        /*
         * Determines the total number of Fluorescence photons detected.
         */
        int NumberFluorescencePhotons(Shower shower);

        /*
         * Determines the total number of Cherenkov photons detected.
         */
        int NumberCherenkovPhotons(Shower shower);

        /*
         * Creates a Cherenkov photon with a randomly-assigned direction (the direction follows a e^-theta/sin(theta)
         * distribution.
         */
        Ray GenerateCherenkovPhoton(Shower shower);

        /*
         * Finds the points where a ray will or has intersected with a sphere centered at the origin. If the ray does
         * not intersect with the sphere, "point" is set to (0, 0, 0) and false is returned. Otherwise, "point" is set
         * to the intersection with the smallest (negative) z-coordinate and "true" is returned.
         */
        bool BackOriginSphereImpact(Ray ray, TVector3* point, double radius);

        /*
         * Determines whether the xy projection of the vector lies within a disk centered at the origin.
         */
        bool WithinXYDisk(TVector3 vec, double radius);

        /*
         * Determines the local fluorescence yield of the shower.
         */
        double FluorescenceYield(Shower shower);

        /*
         * Calculates the effective ionization loss rate for a shower (alpha_eff).
         */
        double IonizationLossRate(Shower shower);

        /*
         * Calculates the atmospheric density at some height above the origin (the location of the detector). The height
         * is NOT relative to sea level.
         */
        double AtmosphereDensity(double height);

        /*
         * Calculates delta = n - 1 for the atmosphere at some height. This assumes that the quantity delta is
         * approximately proportional to the local atmospheric density.
         */
        double AtmosphereDelta(double height);

        /*
         * Calculates the Cherenkov threshold energy at some height. Uses the AtmosphereDelta function.
         */
        double EThresh(double height);

        /*
         * Calculates the critical angle in the expression for the Cherenkov angular distribution.
         */
        double ThetaC(double height);


    public:

        /*
         * Parses the specified file to XML and attemts to extract parameters.
         */
        void ParseFile(boost::property_tree::ptree config);

        PhotonCount SimulateShower(Shower shower);
    };
}

#endif

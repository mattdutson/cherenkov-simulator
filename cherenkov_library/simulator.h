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
#include "common.h"
#include "TF1.h"
#include "Math/Transform3D.h"
#include "TRandom3.h"
#include "TRotation.h"

namespace cherenkov_simulator
{
    class Simulator
    {

    private:

        // The atmospheric scale height
        double H;

        // The density of air at the detector
        double rho_0;

        // 1 - the index of refraction at the detector
        double delta_0;

        // The number of steps to increment the shower
        int n_steps;

        // A general purpose random number generator
        TRandom3 rng;

        double lambda;

        // Various probability distributions
        TF1 energy_distribution;
        TF1 cosine_distribution;
        TF1 impact_distrubition;
        TF1 interact_distribution;

        double stop_radius;

        // A transformation which rotates us from the detector frame to the world frame.
        TRotation rotate_to_world;

        // There are problems with the copy constructor of the FileOptions class (due to the boost library used), so we
        // don't copy it and just use a pointer.
        FileOptions* config;

        void ViewFluorescencePhotons(Shower shower, PhotonCount* photon_count);

        void ViewCherenkovPhotons(Shower shower, Plane ground_plane, PhotonCount* photon_count);

        void SimulateOptics(Ray photon, PhotonCount* photon_count);

        /*
         * Finds the amount that the position of a given shower should be incremented based on the local atmospheric
         * density.
         */
        double FindDistance(Shower shower);

        void AddNoise(PhotonCount* photon_count);

        VoltageSignal VoltageResponse(PhotonCount photon_count);

        TVector3 RandomStopImpact();

        TVector3 MirrorNormal(TVector3 point);

        bool LensImpactPoint(Ray ray, TVector3* point);

        bool MirrorImpactPoint(Ray ray, TVector3* point);

        bool CameraImpactPoint(Ray ray, TVector3* point);

        bool BlockedByCamera(TVector3 start, TVector3 end);

        void DeflectFromLens(Ray* photon);

        bool ReflectFromGround(Ray* photon);

        /*
         * Finds the x and y camera index given a 3-d impact location.
         */
        void ImpactPointToCameraIndex(TVector3 impact, int* x_index, int* y_index);

        /*
         * Finds the direction of a photomultiplier given its x and y indices in the array.
         */
        TVector3 CameraIndexToViewDirection(int x_index, int y_index);

        TVector3 GetViewDirection(TVector3 impact_point);

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
        int NumberFluorescencePhotons(Shower shower, double depth);

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
         * Determines the slant depth between two points (g/cm^2), assuming an exponential atmosphere with parameters specified
         * in the configuration file.
         */
        double SlantDepth(TVector3 point1, TVector3 point2);

        /*
         * Determines the vertical depth beetween two points (g/cm^2), using the same assumptions as SlantDepth.
         */
        double VerticalDepth(TVector3 point1, TVector3 point2);

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
        double Simulator::FluorescenceYield(Shower shower);

        /*
         * Calculates the current age of the shower.
         */
        double ShowerAge(Shower shower);

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

        TVector3 FitSDPlane(PhotonCount data);

    public:

        Simulator(FileOptions* config);

        VoltageSignal SimulateShower(Shower shower);

        Shower ReconstructShower(VoltageSignal dat);

        void EstimateAccuracy();

        Shower GenerateRandomShower();
    };
}

#endif

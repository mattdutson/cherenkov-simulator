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

        // Simulation-life parameters
        double H;
        double rho_0;
        double theta_0;
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

        void ImpactPointToCameraIndex(TVector3 impact, int* x_index, int* y_index);

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
         * Determines the slant depth between two points (g/cm^2), assuming an exponential atmosphere with parameters specified
         * in the configuration file.
         */
        double SlantDepth(TVector3 point1, TVector3 point2);

        /*
         * Determines the vertical depth beetween two points (g/cm^2), using the same assumptions as SlantDepth.
         */
        double VerticalDepth(TVector3 point1, TVector3 point2);

    public:

        Simulator(FileOptions* config);

        VoltageSignal SimulateShower(Shower shower);

        Shower ReconstructShower(VoltageSignal dat);

        void EstimateAccuracy();

        Shower GenerateRandomShower();
    };
}

#endif

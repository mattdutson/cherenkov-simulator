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

namespace cherenkov_simulator
{
    class Simulator
    {

    private:

        // Simulation-life parameters
        double min_energy;
        double max_energy;
        double energy_power;
        double impact_param_min;
        double impact_param_max;
        double H;
        double rho0;
        double detector_elevation;
        int n_steps;

        double lambda;

        // A function of atmosphere
        TF1 energy_distribution;
        TF1 atmosphere;



        Transformation to_detector_frame;

        Transformation to_world_frame;

        FileOptions config;

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

        double FluorescenceFractionCaptured(Shower shower);

        TVector3 MirrorNormal(TVector3 point);

        bool LensImpactPoint(Ray ray, TVector3* point);

        bool MirrorImpactPoint(Ray ray, TVector3* point);

        bool CameraImpactPoint(Ray ray, TVector3* point);

        bool BlockedByCamera(TVector3 start, TVector3 end);

        void DeflectFromLens(Ray* photon);

        bool ReflectFromGround(Ray* photon);

        void ImpactPointToCameraIndex(TVector3 impact, int* x_index, int* y_index);

        TVector3 GetViewDirection(TVector3 impact_point);

        int NumberFluorescencePhotons(Shower shower);

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

        VoltageSignal SimulateShower(Shower shower);

        Shower ReconstructShower(VoltageSignal dat);

        void EstimateAccuracy();

        Shower GenerateRandomShower();
    };
}

#endif

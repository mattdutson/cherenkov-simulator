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

namespace cherenkov_simulator
{
    class Simulator
    {
        
    private:
        
        FileOptions config;
        
        void ViewFluorescencePhotons(Shower shower, PhotonCount* photon_count);
        
        void ViewCherenkovPhotons(Shower shower, Plane ground_plane, PhotonCount* photon_count);
        
        void SimulateOptics(Ray photon, PhotonCount* photon_count);
        
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
        
    public:
        
        VoltageSignal SimulateShower(Shower shower);
        
        Shower ReconstructShower(VoltageSignal dat);
        
        void EstimateAccuracy();
    };
}

#endif

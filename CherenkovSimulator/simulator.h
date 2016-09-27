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
        
        void ViewPoint(Shower shower, RawData* raw_data);
        
        void CherenkovBlast(Shower shower, RawData* raw_data);
        
        SignalData BinRawData(RawData raw_data);
        
        void AddNoise(SignalData* binned);
        
        SignalData CalculateResponse(SignalData count);
        
        TVector3 RandomStopImpact();
        
        double FractionCaptured(Shower shower);
        
        TVector3 OuterLensNormal(TVector3 point);
        
        TVector3 InnerLensNormal(TVector3 point);
        
        TVector3 MirrorNormal(TVector3 point);
        
        bool LensExitPoint(Ray ray, TVector3* point);
        
        bool MirrorImpactPoint(Ray ray, TVector3* point);
        
        bool CameraImpactPoint(Ray ray, TVector3* point);
        
        TVector3 GetViewDirection(TVector3 impact_position);
        
        bool CameraImpact(TVector3 start, TVector3 end);
        
    public:
        
        SignalData SimulateShower(Shower shower);
        
        Shower ReconstructShower(SignalData dat);
        
        void EstimateAccuracy();
    };
}

#endif

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
        
        void ViewPoint(Shower shower, RawData raw_data);
        
        void CherenkovBlast(Shower shower, RawData raw_data);
        
        SignalData BinRawData(RawData raw_data);
        
        void AddNoise(SignalData binned);
        
        SignalData CalculateResponse(SignalData count);
        
    public:
        
        SignalData SimulateShower(Shower shower);
        
        Shower ReconstructShower(SignalData dat);
        
        void EstimateAccuracy();
    };
}

#endif

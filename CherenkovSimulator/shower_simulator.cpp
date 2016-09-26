// shower_simulator.cpp
// cherenkov_simulator
//
// Created by Matthew Dutson on 9/8/16.
//
//

#include "simulator.h"
#include "geometric_objects.h"

namespace cherenkov_simulator
{
    SignalData Simulator::SimulateShower(Shower shower)
    {
        double time_delay = config.Get<double>("time_delay");
        Plane ground_plane = Plane(config.Get<TVector3>("ground_normal"), config.Get<TVector3>("ground_point"));
        int n_steps = (int) (shower.TimeToPlane(ground_plane) / time_delay) + 2;
        
        RawData raw_data = RawData();
        
        for(int i = 0; i < n_steps; i++) {
            ViewPoint(shower, raw_data);
            shower.IncrementPosition(time_delay);
        }
        
        CherenkovBlast(shower, raw_data);
        
        SignalData binned = BinRawData(raw_data);
        
        AddNoise(binned);
        
        SignalData voltage = CalculateResponse(binned);
        
        return voltage;
    }
}

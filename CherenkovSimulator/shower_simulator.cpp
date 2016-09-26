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
        RawData raw_data = RawData();
        
        for(int i = 0; i < n_steps; i++) {
            ViewPoint(shower, raw_data);
            shower.IncrementPosition(timeDelay);
        }
        
        CherenkovBlast(shower, raw_data);
        
        SignalData binned = BinRawData(raw_data);
        
        AddNoise(binned);
        
        SignalData voltage = CalculateResponse(raw_data);
        
        return voltage;
    }
}

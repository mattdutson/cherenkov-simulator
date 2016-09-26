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
        
        ConfigManager config;
        
    public:
        
        SignalData SimulateShower(Shower shower);
    };
}

#endif

// shower_simulator.h
// cherenkov_simulator
//
// Created by Matthew Dutson on 9/8/16.
//
//

#ifndef shower_simulator_h
#define shower_simulator_h

#include "data_containers.h"
#include "geometric_objects.h"

namespace cherenkov_simulator {
   
    namespace shower_simulator {
        
        SignalData SimulateShower(Shower shower);
    }
}

#endif

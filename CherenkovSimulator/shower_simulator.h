/*
 * ShowerSimulator.h
 * CherenkovSimulator
 *
 * Created by Matthew Dutson on 9/8/16.
 *
 *
 */

#ifndef ShowerSimulator_h
#define ShowerSimulator_h

#include "DataContainers.h"
#include "GeometricObjects.h"

namespace shower_simulator {
    
    SignalData SimulateShower(Shower shower);
}

#endif

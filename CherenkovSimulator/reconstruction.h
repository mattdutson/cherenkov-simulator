/*
 * Reconstruction.h
 * CherenkovSimulator
 *
 * Created by Matthew Dutson on 9/8/16.
 *
 *
 */

#ifndef Reconstruction_h
#define Reconstruction_h

#include "DataContainers.h"
#include "GeometricObjects.h"

namespace reconstruction {
    
    Shower ReconstructShower(SignalData dat);
}

#endif

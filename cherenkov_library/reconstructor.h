// reconstructor.h
// cherenkov_library
//
// Created by Matthew Dutson on 1/8/17.
//

#ifndef reconstructor_h
#define reconstructor_h

#include "data_containers.h"
#include "geometric_objects.h"
#include <istream>

namespace cherenkov_simulator
{
    class Reconstructor
    {
    public:

        void ParseFile(std::ifstream config_file);

        TVector3 FitSDPlane(PhotonCount data);

    private:

        // Parameters relating to the position and orientation of the detector relative to its surroundings
        Plane ground_plane;
        TRotation rotate_to_world;

    };
}

#endif

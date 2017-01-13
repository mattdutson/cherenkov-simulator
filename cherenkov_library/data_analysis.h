//
// Created by Matthew Dutson on 1/9/17.
//

#ifndef analysis_h
#define analysis_h

#include "data_containers.h"
#include "geometric_objects.h"

namespace cherenkov_simulator
{
    bool Above(TVector3 reference, TVector3 other);

    std::vector<std::array<double, 2>> CollapseToProfile(PhotonCount data, Plane s_d_plane);
}

#endif
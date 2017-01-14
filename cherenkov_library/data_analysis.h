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

    /*
     * Finds the superimposed time profile of all photomultiplier signals. It is assumed that all signals will have the
     * same time length (if AddNoise wasn't called this won't be the case).
     */
    std::vector<std::array<double, 2>> SuperimposeTimes(PhotonCount data);
}

#endif
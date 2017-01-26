//
// Created by Matthew Dutson on 1/9/17.
//

#ifndef analysis_h
#define analysis_h

#include <TH2C.h>
#include <vector>

#include "data_containers.h"
#include "geometric_objects.h"

namespace cherenkov_library
{
    bool Above(TVector3 reference, TVector3 other);

    std::vector<std::vector<double>> CollapseToProfile(PhotonCount data, Plane s_d_plane);

    /*
     * Finds the superimposed time profile of all photomultiplier signals. It is assumed that all signals will have the
     * same time length (if AddNoise wasn't called this won't be the case).
     */
    std::vector<std::vector<double>> SuperimposeTimes(PhotonCount data);

    /*
     * Creates a 2D histogram with a 1 for valid cells and a 0 for invalid cells.
     */
    TH2C GetValidMap(PhotonCount data);
}

#endif
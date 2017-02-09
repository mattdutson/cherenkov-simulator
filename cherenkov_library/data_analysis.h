//
// Created by Matthew Dutson on 1/9/17.
//

#ifndef analysis_h
#define analysis_h

#include <TH2C.h>
#include <vector>
#include <TGraph.h>

#include "data_containers.h"
#include "geometric_objects.h"

namespace cherenkov_library
{


    void CollapseToProfile(PhotonCount data, Plane s_d_plane, TVector3 shower_axis, std::vector<double>* angles,
                           std::vector<double>* counts);

    /*
     * Finds the superimposed time profile of all photomultiplier signals. It is assumed that all signals will have the
     * same time length (if AddNoise wasn't called this won't be the case).
     */
    void SuperimposeTimes(PhotonCount data, std::vector<double>* times, std::vector<double>* counts);

    /*
     * Creates a 2D histogram with a 1 for valid cells and a 0 for invalid cells.
     */
    TH2C GetValidMap(PhotonCount data);

    /*
     * Makes a TGraph of the collapsed time profile for the simulation.
     */
    TGraph MakeProfileGraph(PhotonCount data);

    /*
     * Makes a 2D histogram with the sum of each bin
     */
    TH2I MakeSumMap(PhotonCount data);
}

#endif
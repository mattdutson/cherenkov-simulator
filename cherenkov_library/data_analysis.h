// data_analysis.h
// cherenkov_library
//
// Created by Matthew Dutson on 1/9/17.
//
// Contains methods which may be useful for analyzing data produced by methods in the Simulator, MonteCarlo, and
// Reconstructor classes. Declared separately from utility.h to prevent circular dependencies.

#ifndef analysis_h
#define analysis_h

#include <TH2C.h>
#include <vector>
#include <TGraph.h>

#include "data_containers.h"
#include "geometric_objects.h"

namespace cherenkov_library
{
    class DataAnalysis
    {

    public:

        /*
         * Creates a 2-d plot of number of photons vs. angle in the shower-detector plane, with one data point for each
         * pixel.
        */
        static void
        CollapseToProfile(PhotonCount data, Plane s_d_plane, TVector3 shower_axis, std::vector<double>* angles,
                          std::vector<double>* counts);

        /*
         * Finds the superimposed time profile of all photomultiplier signals. It is assumed that all signals will have the
         * same time length (if AddNoise wasn't called this won't be the case).
         */
        static void SuperimposeTimes(PhotonCount data, std::vector<double>* times, std::vector<double>* counts);

        /*
         * Creates a 2D histogram with a 1 for valid cells and a 0 for invalid cells.
         */
        static TH2C GetValidMap(PhotonCount data);

        /*
         * Creates a 2D histogram with a 1 for true values and a 0 for false values.
         */
        static TH2C GetBooleanMap(std::vector<std::vector<bool>> valid);

        /*
         * Makes a TGraph of the collapsed time profile for the simulation.
         */
        static TGraph MakeProfileGraph(PhotonCount data);

        /*
         * Makes a 2D histogram with the sum of each bin
         */
        static TH2I MakeSumMap(PhotonCount data);
    };
}

#endif
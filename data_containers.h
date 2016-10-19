// data_containers.h
// cherenkov_simulator
//
// Created by Matthew Dutson on 9/8/16.
//
//

#ifndef data_containers_h
#define data_containers_h

#include "common.h"
#include "TVector3.h"

namespace cherenkov_simulator
{

    class PhotonCount
    {
    private:
        
        std::vector<std::vector<std::vector<int>>> photon_counts;
        
        FileOptions config;
        
    public:
        
        PhotonCount(FileOptions config);
        
        void AddPhoton(double time, int x_index, int y_index);

        void AddPoint(TVector3 direction, double time);
        
    };
    
    class VoltageSignal
    {
        
    };
}

#endif

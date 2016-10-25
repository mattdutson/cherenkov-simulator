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

        class SignalIterator
        {

            friend class PhotonCount;

        private:

            int x_current;

            int y_current;

            std::vector<std::vector<bool>> nonempty_pixels;

            SignalIterator();

            SignalIterator(std::vector<std::vector<bool>> pixel_map);

        public:

            int X();

            int Y();

            bool Next();

            void Reset();

        };
        
        PhotonCount(FileOptions config);
        
        void AddPhoton(double time, int x_index, int y_index);

        void AddPoint(TVector3 direction, double time);

        TVector3 PixelDirection(int x_index, int y_index);

        SignalIterator Iterator();
    };
    
    class VoltageSignal
    {
        
    };
}

#endif

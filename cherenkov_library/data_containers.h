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

        bool ValidPixel(int x_index, int y_index);
        
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

            /*
             * Moves to the next photomultiplier signal. Returns false if the iterator has reached the end of the
             * collection (the last photomultiplier).
             */
            bool Next();

            void Reset();

        };

        /*
         * Creates a 2d indexed collection of vectors (a vector of vectors of vectors). Each vector represents the time
         * series for a particular photomultiplier. The pixel_count parameter defines the width/height of the cluster in
         * number of pixels.
         */
        PhotonCount(int pixel_count);

        void AddPhoton(int time_slot, int x_index, int y_index);

        TVector3 PixelDirection(int x_index, int y_index);

        /*
         * Returns an object for iterating through the pixels.
         */
        SignalIterator Iterator();

        /*
         * Sums all bins in the channel specified by the iterator.
         */
        double SumBins(SignalIterator iter);
    };
    
    class VoltageSignal
    {
        
    };
}

#endif

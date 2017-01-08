// data_containers.h
// cherenkov_simulator
//
// Created by Matthew Dutson on 9/8/16.
//
//

#ifndef data_containers_h
#define data_containers_h

#include "utility.h"
#include "TVector3.h"

namespace cherenkov_simulator
{

    class PhotonCount
    {
    private:

        // TODO: Pass this container information about the relationship between pixel directions and indices.

        std::vector<std::vector<std::vector<int>>> photon_counts;

        // The time at the beginning of the zeroth bin
        double start_time;

        // The size of time bins
        double time_bin;

        // The angle viewed by each photomultiplier
        double pmt_angular_size;

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
        PhotonCount(int n_pmt_across, double start_time, double time_bin, double pmt_angular_size);

        void AddPhoton(double time, TVector3 direction);

        void AddNoise(double noise_rate, SignalIterator current);

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
}

#endif

// DataStructures.h
//
// Author: Matthew Dutson
//
// Contains the definition of the PhotonCount container as well as an iterator to move through it.

#ifndef DATA_STRUCTURES_H
#define DATA_STRUCTURES_H

#include <vector>
#include <TVector3.h>
#include <TRandom3.h>

namespace cherenkov_simulator
{
    /*
     * A class containing a 2d collection of vectors. Each vector represents the photon-time signal for a particular
     * photomultiplier. Also contains basic information about the surroundings which is used to find the direction of
     * each photomultiplier.
     */
    class PhotonCount
    {
    public:

        /*
         * An object used for keeping track of indices in the photon count container and iterating through the non-null
         * photomultipliers.
         */
        class Iterator
        {
        public:

            /*
             * The only constructor. Takes a 2D vector of bools which can be used to determine which pixels are valid. The
             * user is responsible for ensuring that all sub-vectors are non-null.
             */
            Iterator(std::vector<std::vector<bool>> validPixels);

            /*
             * Returns the current x index of the iterator.
             */
            int X() const;

            /*
             * Returns the current y index of the iterator.
             */
            int Y() const;

            /*
             * Moves to the next photomultiplier signal. Returns false if the iterator has reached the end of the
             * collection (the last photomultiplier).
             */
            bool Next();

            /*
             * Returns the iterator to its starting state.
             */
            void Reset();

        private:

            friend class DataStructuresTest;

            // Defines the subset of pixels which the iterator is allowed to stop at
            std::vector<std::vector<bool>> valid;

            // The current indices of the iterator
            int curr_x;
            int curr_y;
        };

        /*
         * A container for PhotonCount parameters to be passed to the constructor.
         */
        struct Params
        {
            int n_pixels;
            double bin_size;
            double angular_size;
            double linear_size;
        };

        /*
         * Creates a 2d indexed collection of vectors (a vector of vectors of vectors). Each vector represents the time
         * series for a particular photomultiplier. Takes the width/height of the array in number of photomultipliers,
         * the starting time for the series, the width of a time bin, and the angular width of a photomultiplier.
         */
        PhotonCount(Params params, double start_time);

        /*
         * Returns the 2D array of valid pixel flags.
         */
        std::vector<std::vector<bool>> GetValid();

        /*
         * Returns the width/height of the 2D array
         */
        int Size();

        /*
         * Returns the total number of time bins in the data structure.
         */
        int NBins();

        /*
         * Returns the size, in seconds, of each time bin.
         */
        double BinSize();

        /*
         * Finds the time corresponding to the specified bin.
         */
        double Time(int bin);

        /*
         * Finds the bin corresponding to some time.
         */
        int Bin(double time);

        /*
         * Determines the direction of the photomultiplier referenced by the iterator.
         */
        TVector3 Direction(const Iterator* iter);

        /*
         * Gets the complete time signal at the current position of the iterator.
         */
        std::vector<int> Signal(const Iterator* iter);

        /*
         * Sums all bins in the channel referenced by the iterator.
         */
        int SumBins(const Iterator* iter);

        /*
         * Finds the average time in the bin referenced by the iterator.
         */
        double AverageTime(const Iterator* iter);

        /*
         * Returns an object for iterating through the pixels.
         */
        Iterator GetIterator();

        /*
         * Increments the photon count at some time for the photomultiplier pointing in the specified direction. If the
         * specified time is before the container's start time, nothing is done.
         */
        void AddPhoton(double time, TVector3 direction, double thinning = 1.0);

        /*
         * Adds some number of noise photons to the channel referenced by the iterator. The noise rate represents the
         * number of photons per second per steradian per square centimeter. These photons are randomly scattered
         * throughout the time bins using the random number generator.
         */
        void AddNoise(double noise_rate, const Iterator* iter, TRandom3* rng);

        /*
         * Subtract the average noise rate from the signal in the pixel specified by the iterator.
         */
        void SubtractNoise(double noise_rate, const Iterator* iter);

        /*
         * Clears any bins in the current pixel which are less than hold_thresh * sigma from zero.
         */
        void ClearNoise(const Iterator* iter, double noise_rate, double hold_thresh);

        /*
         * Returns a vector which contains "true" for each bin in the current pixel which contains more than
         * trigger_thresh * sigma photon counts.
         */
        std::vector<bool> FindTriggers(const Iterator* iter, double noise_rate, double trigger_thresh);

        /*
         * Erases any photon counts which are not in a time bin corresponding to a true value in the input vector.
         */
        void EraseNonTriggered(std::vector<bool> good_bins);

        /*
         * Determines the number of photons per bin observed by a single pixel given the number of photons per second
         * per steradian in a given direction.
         */
        double RealNoiseRate(double universal_rate);

        /*
         * Resizes all channels so they have bins up through the last photon seen and all have the same size.
         */
        void Equalize();

    private:

        friend class DataStructuresTest;

        // The underlying data structure and its validity mask
        std::vector<std::vector<std::vector<int>>> counts;
        std::vector<std::vector<bool>> valid;

        // The number and size of each pixel - cgs, sr
        int n_pixels;
        double angular_size;
        double linear_size;

        // The time at the beginning of the zeroth bin - cgs
        double bin_size;
        double start_time;
        double last_time;

        // Keeps track of whether we need to call Equalize
        bool equalized;

        /*
         * Determines whether the pixel at the specified indices lies within the central circle.
         */
        bool ValidPixel(int x_index, int y_index);

        /*
         * A private method which is functionally equivalent to Direction(const Iterator*).
         */
        TVector3 Direction(int x_index, int y_index);

        /*
         * If the vector at index (x, y) is smaller than the specified size, it is expanded. This method does NOT
         * check that the indices are valid.
         */
        void ExpandVector(int x_index, int y_index, int size);
    };
}

#endif

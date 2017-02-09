// data_containers.h
// cherenkov_library
//
// Created by Matthew Dutson on 9/8/16.
//
// Contains the definition of the photon count container as well as an iterator to move through it.

#ifndef data_containers_h
#define data_containers_h

#include <TVector3.h>
#include <TRandom3.h>
#include <vector>

namespace cherenkov_library
{
    /*
     * An object used for keeping track of indices in the photon count container and iterating through the non-null
     * photomultipliers.
     */
    class SignalIterator
    {

    public:

        /*
         * The only constructor. Takes a 2D vector of bools which can be used to determine which pixels are valid. The
         * user is responsible for ensuring that all sub-vectors are non-null.
         */
        SignalIterator(std::vector<std::vector<bool>> validPixels);

        /*
         * Returns the current x index of the iterator.
         */
        int X();

        /*
         * Returns the current y index of the iterator.
         */
        int Y();

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

        // A 2D vector of bools which can be used to determine which pixels are valid.
        std::vector<std::vector<bool>> valid_pixels;

        // The current x index of the iterator
        int x_current;

        // The current y index of the iterator
        int y_current;
    };

    /*
     * A class containing a 2d collection of vectors. Each vector represents the photon-time signal for a particular
     * photomultiplier. Also contains basic information about the surroundings which is used to find the direction of
     * each photomultiplier.
     */
    class PhotonCount
    {

    public:

        /*
         * Creates a 2d indexed collection of vectors (a vector of vectors of vectors). Each vector represents the time
         * series for a particular photomultiplier. Takes the width/height of the array in number of photomultipliers,
         * the starting time for the series, the width of a time bin, and the angular width of a photomultiplier.
         */
        PhotonCount(int n_pmt_across, double start_time, double time_bin, double pmt_angular_size,
                    double pmt_linear_size);

        /*
         * Returns the width/height of the 2D array
         */
        int Size();

        /*
         * Returns the total number of bins in the data structure.
         */
        int NBins();

        /*
         * Finds the time corresponding to the specified bin.
         */
        double Time(int bin);

        /*
         * Finds the bin corresponding to some time.
         */
        int Bin(double time);

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
        void AddNoise(double noise_rate, SignalIterator current, TRandom3 rng);

        /*
         * Clears any bins in the current pixel which are less than hold_thresh * sigma from the mean.
         */
        void ClearNoise(SignalIterator current, double noise_rate, double hold_thresh);

        /*
         * Returns a vector which contains "true" for each bin in the current pixel which contains more than
         * trigger_thresh * sigma photon counts.
         */
        std::vector<bool> FindTriggers(SignalIterator current, double noise_rate, double trigger_thresh);

        /*
         * Erases any photon counts which are not in a time bin corresponding to a true value in the input vector.
         */
        void EraseNonTriggered(std::vector<bool> good_bins);

        /*
         * Determines the direction of the photomultiplier referenced by the iterator.
         */
        TVector3 Direction(SignalIterator current);

        /*
         * Gets the complete time signal at the current position of the iterator.
         */
        std::vector<int> Signal(SignalIterator current);

        /*
         * Returns the 2D array of valid pixel flags.
         */
        std::vector<std::vector<bool>> GetValid();

        /*
         * Returns an object for iterating through the pixels.
         */
        SignalIterator Iterator();

        /*
         * Sums all bins in the channel referenced by the iterator.
         */
        int SumBins(SignalIterator iter);

        /*
         * Finds the average time in the bin referenced by the iterator.
         */
        double AverageTime(SignalIterator iter);

    private:

        /*
         * Determines whether the pixel at the specified indices lies within the central circle.
         */
        bool ValidPixel(int x_index, int y_index);

        /*
         * A private method which determines the direction of the photomultiplier at the specified indices.
         */
        TVector3 Direction(int x_index, int y_index);

        /*
         * If the vector at indices (x, y) is smaller than the specified size, it is expanded. This method does NOT
         * check that the indices are valid.
         */
        void ExpandVector(int x_index, int y_index, int min_size);

        /*
         * Determines the number of photons per bin observed by a single pixel given the number of photons per square
         * centimeter per steradian in a given direction.
         */
        double RealNoiseRate(double universal_rate);

        // The underlying data structure to store photon counts
        std::vector<std::vector<std::vector<int>>> photon_counts;

        // A 2D vector of bools which can be used to determine which pixels are valid.
        std::vector<std::vector<bool>> valid_pixels;

        // The number of pixels across
        int n_pixels;

        // The time at the beginning of the zeroth bin
        double start_time;

        // The last real signal (added by AddPhoton) seen by the array
        double last_time;

        // The size of time bins
        double bin_size;

        // The angle viewed by each photomultiplier
        double pmt_angular_size;

        // The side length of each photomultiplier
        double pmt_linear_size;
    };
}

#endif

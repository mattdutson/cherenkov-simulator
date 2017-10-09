// DataStructures.h
//
// Author: Matthew Dutson
//
// Contains the definition of the PhotonCount container as well as an iterator to move through it.

#ifndef DATA_STRUCTURES_H
#define DATA_STRUCTURES_H

#include <vector>
#include <TVector3.h>
#include <TF1.h>

#include "Utility.h"

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
            explicit Iterator(Bool2D validPixels);

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
             * collection (the last photomultiplier). Steps through y and then through x.
             */
            bool Next();

            /*
             * Returns the iterator to its starting state.
             */
            void Reset();

        private:

            friend class DataStructuresTest;

            // Defines the subset of pixels which the iterator is allowed to stop at
            Bool2D valid;

            // The current indices of the iterator
            int curr_x;
            int curr_y;
        };

        /*
         * A container for PhotonCount parameters to be passed to the constructor.
         */
        struct Params
        {
            size_t n_pixels;
            size_t max_bytes;
            double bin_size;
            double angular_size;
            double linear_size;
        };

        /*
         * Default constructor
         */
        PhotonCount();

        /*
         * Creates a 2d indexed collection of vectors (a vector of vectors of vectors). Each vector represents the time
         * series for a particular photomultiplier. Takes the width/height of the array in number of photomultipliers,
         * the starting time for the series, the width of a time bin, and the angular width of a photomultiplier. The
         * width/height of the array can be zero (this is what the default constructor does), although this probably
         * isn't very useful. The bin size CANNOT be zero, as this results in divisions by zero.
         */
        PhotonCount(Params params, double min_time, double max_time);

        /*
         * Returns the 2D array of valid pixel flags.
         */
        Bool2D GetValid() const;

        /*
         * Returns the width/height of the 2D array
         */
        size_t Size() const;

        /*
         * Returns the total number of time bins in the data structure.
         */
        size_t NBins() const;

        /*
         * Returns true if no photons were added.
         */
        bool Empty() const;

        /*
         * Finds the time corresponding to the specified bin. Throws an exception if the bin is out of range.
         */
        double Time(int bin) const;

        /*
         * Finds the bin corresponding to some time. Throws an exception if the time is out of range.
         */
        size_t Bin(double time) const;

        /*
         * The angular size of the detector, from axis to the outside of the detector surface.
         */
        double DetectorAxisAngle() const;

        /*
         * Determines the direction of the photomultiplier referenced by the iterator. If the iterator does not have the
         * same size as this object's validity mask, throws an invalid_argument exception.
         */
        TVector3 Direction(const Iterator& iter) const;

        /*
         * Gets the complete time signal at the current position of the iterator.
         */
        Short1D Signal(const Iterator& iter) const;

        /*
         * Sums all bins in the channel referenced by the iterator.
         */
        int SumBins(const Iterator& iter) const;

        /*
         * Sums all bins in the channel which correspond to a "true" value in the filter.
         */
        int SumBinsFiltered(const Iterator& iter, const Bool3D* filter) const;

        /*
         * Finds the average time in the pixel referenced by the iterator. Throws a domain_error if the channel
         * referenced by the iterator is empty (this results in division by zero).
         */
        double AverageTime(const Iterator& iter) const;

        /*
         * Finds the standard deviation of the times in the pixel referenced by the iterator, applying Sheppard's
         * correction because the data is binned. Throws a domain_error if the channel referenced by the iterator is
         * empty (this results in division by zero).
         */
        double TimeError(const Iterator& iter) const;

        /*
         * Returns an object for iterating through the pixels.
         */
        Iterator GetIterator() const;

        /*
         * Returns a 3D matrix of false values, with the same dimensions as the underlying vector of the PhotonCount
         * class.
         */
        Bool3D GetFalseMatrix() const;

        /*
         * Increments the photon count at some time for the photomultiplier pointing in the specified direction. If the
         * specified time is before the container's start time, nothing is done.
         */
        void AddPhoton(double time, TVector3 direction, int thinning);

        /*
         * Adds some number of noise photons to the channel referenced by the iterator. The noise rate represents the
         * number of photons per second per steradian per square centimeter. These photons are randomly scattered
         * throughout the time bins using the random number generator. Adding noise will set "empty" to false,
         * regardless of how much (if any) noise is added.
         */
        void AddNoise(double noise_rate, const Iterator& iter);

        /*
         * Subtract the average noise rate from the signal in the pixel specified by the iterator.
         */
        void Subtract(const Iterator& iter, double rate);

        /*
         * Erases any photon counts which do not correspond to a true value in the 3D input vector.
         */
        void Subset(const Bool3D& good_pixels);

        /*
         * Returns a vector which contains "true" for each bin in the current pixel which contains more than
         * trigger_thresh * sigma photon counts.
         */
        Bool1D AboveThreshold(const Iterator& iter, int threshold) const;

        /*
         * Determines the appropriate threshold given the global noise rate (Poisson distributed), and the number of
         * standard deviations above the mean where the threshold should lie. Any signals at or above this threshold are
         * considered "good".
         */
        int FindThreshold(double global_rate, double sigma) const;

        /*
         * Resizes all channels so they have bins up through the last photon seen and all have the same size.
         */
        void Trim();

    private:

        friend class DataStructuresTest;

        // The underlying data structure and its validity mask
        // TODO: The three-dimensional array wastes some space (fitting a circle into a square) - it would be nice to convert x-y coordinates to a single array index
        Short3D counts;
        Short2D sums;
        Bool2D valid;

        // The number and size of each pixel - cgs, sr
        size_t n_pixels;
        double angular_size;
        double linear_size;

        // The time at the beginning of the zeroth bin - cgs
        double bin_size;
        double min_time;
        double max_time;
        double first_time;
        double last_time;

        // Keeps track of whether we need to call Trim
        bool empty;
        bool trimmed;

        /*
         * Modifies some (x, y, t) cell by the specified amount, taking an Iterator instead of an x, y position.
         */
        void IncrementCell(int inc, const Iterator& iter, size_t t);

        /*
         * Modifies some (x, y, t) by the specified amount.
         */
        void IncrementCell(int inc, size_t x_index, size_t y_index, size_t t);

        /*
         * Determines whether the pixel at the specified indices lies within the central circle.
         */
        bool ValidPixel(int x_index, int y_index) const;

        /*
         * A private method which is functionally equivalent to Direction(const Iterator*).
         */
        TVector3 Direction(int x_index, int y_index) const;

        /*
         * Determines the number of photons per bin observed by a single pixel given the number of photons per second
         * per steradian in a given direction.
         */
        double RealNoiseRate(double universal_rate) const;

        /*
         * Approximates the sum of the specified Poisson distribution from some minimum value to infinity. The summation
         * stops when the value of the distribution is less than 10^-6 of the value at the starting point.
         */
        static double PoissonSum(double mean, int min);

        /*
         * Calculates the value of the Poisson distribution with specified mean at the specified value.
         */
        static double Poisson(double mean, int x);
    };
}

#endif

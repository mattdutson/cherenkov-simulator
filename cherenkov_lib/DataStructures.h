// DataStructures.h
//
// Author: Matthew Dutson
//
// Defines PhotonCount and related classes

#ifndef DATA_STRUCTURES_H
#define DATA_STRUCTURES_H

#include <vector>
#include <TVector3.h>

#include "Utility.h"

namespace cherenkov_simulator
{
    /*
     * A class containing a 2D collection of vectors. Each vector is a histogram of photon arrival times for a
     * particular photomultiplier. Also contains basic information about the detector which is used to find the
     * direction of each photomultiplier.
     */
    class PhotonCount
    {
    public:

        /*
         * An object used for iterating through the valid circular subset of pixels (our 2D data structure places a
         * circle inside a square, so some vectors don't correspond to valid pixels).
         */
        class Iterator
        {
        public:

            /*
             * The only constructor. Takes a 2D vector of booleans with true values for valid pixels. The user is
             * responsible for ensuring that all sub-vectors are non-null.
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
             * Moves to the next valid pixel. Returns false if the iterator has reached the end of the collection.
             * Steps through y first and then through x.
             */
            bool Next();

            /*
             * Reverts the iterator to its starting state.
             */
            void Reset();

        private:

            friend class DataStructuresTest;

            // Contains true values for valid pixels.
            Bool2D valid;

            int curr_x;
            int curr_y;
        };

        /*
         * A container for PhotonCount constructor parameters.
         */
        struct Params
        {
            size_t n_pixels;
            size_t max_byte;
            double bin_size;
            double ang_size;
            double lin_size;
        };

        /*
         * The default constructor. Objects constructed with this should only be used as placeholders.
         */
        PhotonCount();

        /*
         * The main constructor. Takes the size of the array, the maximum amount of memory available, the time bin size,
         * and the size of each individual pixel. Also takes upper and lower limits on the arrival times of photons.
         * Throws an invalid_argument exception if any parameters are out of range.
         */
        PhotonCount(Params params, double min_time, double max_time);

        /*
         * Returns a 2D vector of booleans with true values for valid pixels.
         */
        Bool2D GetValid() const;

        /*
         * Returns the diameter of the pixel array in number of pixels. This is also the size of the underlying 2D array
         * of vectors.
         */
        size_t Size() const;

        /*
         * Returns the number of bins in each of the time series vectors. This can change if Trim() is called.
         */
        size_t NBins() const;

        /*
         * Returns true if no photons have been added. Emptiness will not change in the AddNoise method if all the
         * random noise values are zero.
         */
        bool Empty() const;

        /*
         * Finds the time corresponding to the specified bin. Throws an out_of_range exception if the bin is not in the
         * range [0, NBins()-1].
         */
        double Time(int bin) const;

        /*
         * Finds the bin corresponding to the specified time. Throws an out_of_range exception if the time is not in the
         * range [min_time, max_time], where min_time and max_time were the values passed to the constructor. The valid
         * range may change if Trim() is called.
         */
        size_t Bin(double time) const;

        /*
         * Returns the angle from the axis to the outside of the field of view.
         */
        double DetectorAxisAngle() const;

        /*
         * Determines the direction seen by the pixel at the current location of the iterator.
         */
        TVector3 Direction(const Iterator& iter) const;

        /*
         * Returns the 1D histogram of photon arrival times at the current location of the iterator.
         */
        Short1D Signal(const Iterator& iter) const;

        /*
         * Sums the 1D vector at the current location of the iterator.
         */
        int SumBins(const Iterator& iter) const;

        /*
         * Sums the bins of the 1D vector at the current location of the iterator which correspond to "true" values in
         * the filter.
         */
        int SumBinsFiltered(const Iterator& iter, const Bool3D& filter) const;

        /*
         * Finds the average time in the pixel referenced by the iterator. Throws a domain_error exception if
         * SumBins(iter) == 0 (this results in division by zero).
         */
        double AverageTime(const Iterator& iter) const;

        /*
         * Finds the standard deviation of the mean calculated by AverageTime(). A Sheppard correction is added due to
         * the binned nature of the data. Throws a domain_error exception if SumBins(iter) == 0 (this results in
         * division by zero).
         */
        double TimeError(const Iterator& iter) const;

        /*
         * Returns a PhotonCount::Iterator for this PhotonCount object.
         */
        Iterator GetIterator() const;

        /*
         * Returns a 3D vector of false values with the same dimensions as the PhotonCount underlying 3D vector.
         */
        Bool3D GetFalseMatrix() const;

        /*
         * Increments a bin of the photon count histogram of the pixel at the specified position. Nothing is done
         * if the time is outside [min_time, max_time] or the position is outside the disk-shaped pixel array. An
         * invalid_argument exception is thrown if the position vector is zero. Note that the position specified for
         * this method is not the same as the return from Direction(). direction = -position.Unit().
         */
        void AddPhoton(double time, TVector3 position, int thinning);

        /*
         * Adds background noise to the time series at the specified position. A random Poisson value is generated for
         * each bin at this position. The input noise rate is the number per second per steradian. This is converted to
         * an average rate for a single time bin of a pixel.
         */
        void AddNoise(double noise_rate, const Iterator& iter);

        /*
         * Subtract the average noise rate from the signal in the pixel specified by the iterator. Like AddNoise(), the
         * input rate is the number per second per steradian, which is converted to the appropriate units.
         * The most probable value for a Poisson distribution with mean < 1 is zero, so any non-integer rate is floored.
         */
        void Subtract(double noise_rate, const Iterator& iter);

        /*
         * Returns a vector which contains "true" for each bin in the specified pixel which contains more than a certain
         * number of photons.
         */
        Bool1D AboveThreshold(const Iterator& iter, int threshold) const;

        /*
         * Determines the appropriate threshold given the noise rate (in number per second per sr per square cm) and the
         * number of standard deviations above the mean where the threshold should be set. This is done by upping the
         * threshold until the sum of all Poisson probabilities above the threshold is less than the integral of a
         * Gaussian above the sigma multiple.
         */
        int FindThreshold(double noise_rate, double sigma) const;

        /*
         * Zeroes any photon counts which do not correspond to a true value in the 3D input vector.
         */
        void Subset(const Bool3D& good_pixels);

        /*
         * Resizes all 1D count vectors to remove any leading or trailing segments which are empty in all pixels.
         */
        void Trim();

    private:

        friend class DataStructuresTest;

        Short3D counts;
        Short2D sums;
        Bool2D valid;

        // The number and size of pixels - cgs, sr.
        size_t n_pixels;
        double ang_size;
        double lin_size;

        // Various properties of the time series - cgs.
        double bin_size;
        double min_time;
        double max_time;
        double frst_time;
        double last_time;

        bool empty;
        bool trimd;

        /*
         * A wrapper to the other IncrementCell method which takes an Iterator instead of an (x, y) coordinate.
         */
        void IncrementCell(int inc, const Iterator& iter, size_t t);

        /*
         * Modifies some (x, y, t) bin by the specified amount. Updates the sum of that pixel and changes the empty flag
         * if needed.
         */
        void IncrementCell(int inc, size_t x_index, size_t y_index, size_t t);

        /*
         * Determines whether the pixel at the specified indices lies within the central circle.
         */
        bool IsValid(int x_index, int y_index) const;

        /*
         * A private method which is functionally equivalent to Direction(const Iterator*).
         */
        TVector3 Direction(int x_index, int y_index) const;

        /*
         * Determines the average number of noise photons per bin in a single pixel from the noise rate in number per
         * second per steradian.
         */
        double RealNoiseRate(double noise_rate) const;

        /*
         * Computes the sum of a Poisson distribution with specified mean from the min value (inclusive) to infinity.
         */
        static double PoissonSum(double mean, int min);

        /*
         * Calculates a particular value of a Poisson distribution with specified mean.
         */
        static double Poisson(double mean, int x);
    };
}

#endif

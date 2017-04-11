// DataStructures.cpp
//
// Author: Matthew Dutson
//
// Implementation of DataStructures.h

#include <TMath.h>

#include "DataStructures.h"

using namespace TMath;

using std::vector;

namespace cherenkov_simulator
{
    PhotonCount::Iterator::Iterator(Bool2D validPixels)
    {
        this->valid = validPixels;
        Reset();
    }

    size_t PhotonCount::Iterator::X() const
    {
        if (curr_y < 0) throw std::runtime_error("Invalid position");
        return (size_t) curr_x;
    }

    size_t PhotonCount::Iterator::Y() const
    {
        if (curr_y < 0) throw std::runtime_error("Invalid position");
        return (size_t) curr_y;
    }

    bool PhotonCount::Iterator::Next()
    {
        // If the outer vector is empty, then there will never be anything to iterate through.
        if (valid.size() == 0) return false;

        // Loop until we find a valid pixel.
        bool found = false;
        while (!found)
        {
            // If we're at the end of the current column, try to move to the next column. Otherwise, increment y.
            if (curr_y == valid[curr_x].size() - 1)
            {
                if (curr_x == valid.size() - 1) return false;
                curr_x++;
                curr_y = -1;
            }
            else
            {
                curr_y++;
            }
            if (curr_y > -1) found = valid[curr_x][curr_y];
        }
        return true;
    }

    void PhotonCount::Iterator::Reset()
    {
        curr_x = 0;
        curr_y = -1;
    }

    PhotonCount::PhotonCount(Params params, double min_time, double max_time)
    {
        // Copy parameters
        n_pixels = params.n_pixels;
        bin_size = params.bin_size;
        angular_size = params.angular_size;
        linear_size = params.linear_size;
        this->min_time = min_time;
        this->max_time = max_time;

        // We haven't seen any non-noise photons yet, so the last time is the start time and channels are trimmed
        first_time = max_time;
        last_time = min_time;
        trimmed = false;
        empty = true;

        // Initialize the photon count and valid pixel structures.
        counts = Int3D(n_pixels, Int2D(n_pixels, Int1D()));
        valid = Bool2D(n_pixels, Bool1D(n_pixels, false));
        for (size_t i = 0; i < n_pixels; i++)
        {
            for (size_t j = 0; j < n_pixels; j++)
            {
                valid[i][j] = ValidPixel(i, j);
                if (valid[i][j]) counts[i][j] = Int1D(NBins(), 0);
            }
        }
    }

    Bool2D PhotonCount::GetValid() const
    {
        return valid;
    }

    size_t PhotonCount::Size() const
    {
        return n_pixels;
    }

    size_t PhotonCount::NBins() const
    {
        return min_time == max_time ? 0 : Bin(max_time) + 1;
    }

    bool PhotonCount::Empty() const
    {
        return empty;
    }

    double PhotonCount::Time(int bin) const
    {
        return bin * bin_size + min_time;
    }

    size_t PhotonCount::Bin(double time) const
    {
        return (size_t) Floor((time - min_time) / bin_size);
    }

    double PhotonCount::DetectorAxisAngle() const
    {
        return ((double) n_pixels / 2) * angular_size;
    }

    TVector3 PhotonCount::Direction(const Iterator& iter) const
    {
        return Direction(iter.X(), iter.Y());
    }

    Int1D PhotonCount::Signal(const Iterator& iter) const
    {
        return counts[iter.X()][iter.Y()];
    }

    int PhotonCount::SumBins(const Iterator& iter, const Bool1D* mask) const
    {
        int sum = 0;
        for (size_t i = 0; i < counts[iter.X()][iter.Y()].size(); i++)
        {
            int count = counts[iter.X()][iter.Y()][i];
            if (mask != nullptr) sum += mask->at(i) ? count : 0;
            else sum += count;
        }
        return sum;
    }

    double PhotonCount::AverageTime(const Iterator& iter) const
    {
        int sum = SumBins(iter);
        double average = 0;
        for (int i = 0; i < counts[iter.X()][iter.Y()].size(); i++)
        {
            average += Time(i) * counts[iter.X()][iter.Y()][i] / (double) sum;
        }
        return average;
    }

    double PhotonCount::TimeError(const Iterator& iter) const
    {
        double mean = AverageTime(iter);
        double variance = 0;
        for (int i = 0; i < NBins(); i++)
        {
            variance += counts[iter.X()][iter.Y()][i] * Sq(Time(i) - mean);
        }
        double count = SumBins(iter);
        variance /= count;
        variance += Sq(bin_size) / 12.0;
        return Sqrt(variance / count);
    }

    PhotonCount::Iterator PhotonCount::GetIterator() const
    {
        return Iterator(valid);
    }

    Bool3D PhotonCount::GetFalseMatrix() const
    {
        size_t size = Size();
        return Bool3D(size, Bool2D(size, Bool1D(NBins(), false)));
    }

    void PhotonCount::AddPhoton(double time, TVector3 direction, double thinning)
    {
        if (time < min_time || time > max_time) return;

        // Determine the indices of the pixel based on its orientation
        double elevation = ATan(direction.Y() / direction.Z());
        double azimuth = ATan(direction.X() / direction.Z());
        size_t y_index = (size_t) (Floor(elevation / angular_size) + n_pixels / 2);
        size_t x_index = (size_t) (Floor(azimuth / (angular_size * Cos(elevation))) + n_pixels / 2);

        // If the location is valid, add the pixel to the underlying data structure
        if (ValidPixel(x_index, y_index))
        {
            size_t time_slot = Bin(time);
            counts[x_index][y_index][time_slot] += thinning;
            if (time > last_time) last_time = time;
            if (time < first_time) first_time = time;
            trimmed = false;
            empty = false;
        }
    }

    void PhotonCount::AddNoise(double noise_rate, const Iterator& iter, TRandom3& rng)
    {
        // Ensure all channels have the same length.
        Trim();

        // Randomly determine the number of photons in each bin according to the Poisson distribution.
        double expected_photons = RealNoiseRate(noise_rate);
        for (size_t i = 0; i < NBins(); i++)
        {
            counts[iter.X()][iter.Y()][i] += rng.Poisson(expected_photons);
        }

        empty = false;
    }

    void PhotonCount::Subtract(const Iterator& iter, double rate)
    {
        // We floor the real noise rate because the most probably value for a Poisson distribution with mean < 1 is 0.
        // With the current configuration, mean < 1 seems like a reasonable assumption.
        int expected_photons = (int) RealNoiseRate(rate);
        for (size_t i = 0; i < counts[iter.X()][iter.Y()].size(); i++)
        {
            counts[iter.X()][iter.Y()][i] -= expected_photons;
        }
    }

    Bool1D PhotonCount::AboveThreshold(const Iterator& iter, int threshold) const
    {
        Int1D data = counts[iter.X()][iter.Y()];
        Bool1D triggers = Bool1D(data.size());
        for (size_t i = 0; i < data.size(); i++) triggers[i] = data[i] > threshold;
        return triggers;
    }

    int PhotonCount::FindThreshold(double global_rate, double sigma) const
    {
        double max_prob = Sqrt(Pi() / 2.0) * Erfc(sigma / Sqrt(2));
        double mean = RealNoiseRate(global_rate);
        int thresh = (int) Floor(sigma * Sqrt(mean));
        while(PoissonSum(mean, thresh) > max_prob) thresh++;
        return thresh;
    }

    void PhotonCount::Subset(const Bool3D& good_bins)
    {
        for (size_t i = 0; i < counts.size(); i++)
        {
            for (size_t j = 0; j < counts[i].size(); j++)
            {
                for (size_t t = 0; t < counts[i][j].size(); t++)
                {
                    if (!good_bins[i][j][t]) counts[i][j][t] = 0;
                }
            }
        }
    }

    void PhotonCount::Trim()
    {
        if (trimmed || empty) return;
        Iterator iter = GetIterator();
        while (iter.Next()) {
            Int1D::iterator begin = counts[iter.X()][iter.Y()].begin();
            Int1D::iterator end = counts[iter.X()][iter.Y()].end();
            counts[iter.X()][iter.Y()].erase(begin + Bin(last_time) + 1, end);
            counts[iter.X()][iter.Y()].erase(begin, begin + Bin(first_time));
        }
        min_time = min_time + Floor((first_time - min_time) / bin_size) * bin_size;
        max_time = last_time;
        trimmed = true;
    }

    bool PhotonCount::ValidPixel(size_t x_index, size_t y_index) const
    {
        size_t n_half = n_pixels / 2;
        double arc_length = n_half * linear_size;
        double angle = n_half * angular_size;
        double radius = arc_length / angle;
        double max_dev = radius * Sin(angle);
        TVector3 direction = Direction(x_index, y_index) * radius;
        return Utility::WithinXYDisk(direction, max_dev);
    }

    TVector3 PhotonCount::Direction(size_t x_index, size_t y_index) const
    {
        // Reverse the operations of AddPhoton. Add 0.5 to account for the Floor() operation.
        double pixels_up = (y_index - n_pixels / 2.0 + 0.5);
        double elevation = pixels_up * angular_size;
        double pixels_horiz = (x_index - n_pixels / 2.0 + 0.5);
        double azimuth = pixels_horiz * angular_size * Cos(elevation);

        // A positive azimuth should correspond to a positive x component.
        return TVector3(Cos(elevation) * Sin(azimuth), Sin(elevation), Cos(elevation) * Cos(azimuth));
    }

    double PhotonCount::RealNoiseRate(double universal_rate) const
    {
        // Calculate the number of noise photons per second from the number of photons per second per steradian.
        double solid_angle = Sq(angular_size);
        double pixel_rate = universal_rate * solid_angle;

        // Determine the mean number of Poisson events in a single time bin.
        return bin_size * pixel_rate;
    }

    double PhotonCount::PoissonSum(double mean, int min)
    {
        double sum = 0.0;
        for (int i = min; Poisson(mean, i) > 10e-6 * Poisson(mean, min); i++)
        {
            sum += Poisson(mean, i);
        }
        return sum;
    }

    double PhotonCount::Poisson(double mean, int x)
    {
        return Exp(-mean) * Power(mean, x) / Factorial(x);
    }
}

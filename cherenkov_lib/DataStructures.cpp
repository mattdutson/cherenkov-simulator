// DataStructures.cpp
//
// Author: Matthew Dutson
//
// Implementation of DataStructures.h

#include <TMath.h>

#include "DataStructures.h"
#include "Utility.h"

using namespace TMath;

using std::vector;

namespace cherenkov_simulator
{
    PhotonCount::Iterator::Iterator(std::vector<std::vector<bool>> validPixels)
    {
        this->valid = validPixels;
        Reset();
    }

    int PhotonCount::Iterator::X() const
    {
        return curr_x;
    }

    int PhotonCount::Iterator::Y() const
    {
        return curr_y;
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

    PhotonCount::PhotonCount(Params params, double start_time)
    {
        // Copy parameters
        this->n_pixels = params.n_pixels;
        this->bin_size = params.bin_size;
        this->angular_size = params.angular_size;
        this->linear_size = params.linear_size;
        this->start_time = start_time;

        // We haven't seen any non-noise photons yet, so the last time is the start time and channels are equalized
        last_time = start_time;
        equalized = true;

        // Initialize the photon count and valid pixel structures.
        counts = vector<vector<vector<int>>>();
        valid = vector<vector<bool>>();
        for (int i = 0; i < n_pixels; i++)
        {
            counts.push_back(vector<vector<int>>());
            valid.push_back(vector<bool>());
            for (int j = 0; j < n_pixels; j++)
            {
                counts[i].push_back(vector<int>());
                valid[i].push_back(ValidPixel(i, j));
            }
        }

        // Initialize the gaussian used to calculate threhsolds.
        gauss = TF1("gauss", "e^(-0.5 * x^2) / sqrt(2 * pi)", -Infinity(), Infinity());
    }

    vector<vector<bool>> PhotonCount::GetValid() const
    {
        return valid;
    }

    int PhotonCount::Size() const
    {
        return n_pixels;
    }

    int PhotonCount::NBins() const
    {
        return start_time == last_time ? 0 : Bin(last_time) + 1;
    }

    bool PhotonCount::Empty() const
    {
        return NBins() == 0;
    }

    double PhotonCount::BinSize() const
    {
        return bin_size;
    }

    double PhotonCount::Time(int bin) const
    {
        return bin * bin_size + start_time;
    }

    int PhotonCount::Bin(double time) const
    {
        return (int) Floor((time - start_time) / bin_size);
    }

    double PhotonCount::DetectorAxisAngle() const
    {
        return ((double) n_pixels / 2) * angular_size;
    }

    TVector3 PhotonCount::Direction(const Iterator* iter) const
    {
        return Direction(iter->X(), iter->Y());
    }

    vector<int> PhotonCount::Signal(const Iterator* iter) const
    {
        return counts[iter->X()][iter->Y()];
    }

    int PhotonCount::SumBins(const Iterator* iter, const std::vector<bool>* mask) const
    {
        int sum = 0;
        for (int i = 0; i < counts[iter->X()][iter->Y()].size(); i++)
        {
            int count = counts[iter->X()][iter->Y()][i];
            if (mask != nullptr) sum += mask->at(i) ? count : 0;
            else sum += count;
        }
        return sum;
    }

    double PhotonCount::AverageTime(const Iterator* iter) const
    {
        int sum = SumBins(iter);
        double average = 0;
        for (int i = 0; i < counts[iter->X()][iter->Y()].size(); i++)
        {
            average += Time(i) * counts[iter->X()][iter->Y()][i] / (double) sum;
        }
        return average;
    }

    double PhotonCount::TimeError(const Iterator* iter) const
    {
        double mean = AverageTime(iter);
        double variance = 0;
        for (int i = 0; i < NBins(); i++)
        {
            variance += counts[iter->X()][iter->Y()][i]* Sq(Time(i) - mean);
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

    vector<vector<vector<bool>>> PhotonCount::GetFalseMatrix() const
    {
        int size = Size();
        return vector<vector<vector<bool>>>(size, vector<vector<bool>>(size, vector<bool>(NBins(), false)));
    }

    void PhotonCount::AddPhoton(double time, TVector3 direction, double thinning)
    {
        if (time < start_time) return;

        // Determine the indices of the pixel based on its orientation
        double elevation = ATan(direction.Y() / direction.Z());
        double azimuth = ATan(direction.X() / direction.Z());
        int y_index = (int) (Floor(elevation / angular_size) + n_pixels / 2);
        int x_index = (int) (Floor(azimuth / (angular_size * Cos(elevation))) + n_pixels / 2);

        // If the location is valid, add the pixel to the underlying data structure
        if (ValidPixel(x_index, y_index))
        {
            int time_slot = Bin(time);
            ExpandVector(x_index, y_index, time_slot + 1);
            counts[x_index][y_index][time_slot] += thinning;
            if (time > last_time) last_time = time;
            equalized = false;
        }
    }

    void PhotonCount::AddNoise(double noise_rate, const Iterator* iter, TRandom3* rng)
    {
        // Ensure all channels have the same length.
        Equalize();

        // Randomly determine the number of photons in each bin according to the Poisson distribution.
        double expected_photons = RealNoiseRate(noise_rate);
        for (int i = 0; i < NBins(); i++)
        {
            counts[iter->X()][iter->Y()][i] += rng->Poisson(expected_photons);
        }
    }

    void PhotonCount::Subtract(const Iterator *iter, double rate)
    {
        // We floor the real noise rate because the most probably value for a Poisson distribution with mean < 1 is 0.
        // With the current configuration, mean < 1 seems like a reasonable assumption.
        int expected_photons = (int) RealNoiseRate(rate);
        for (int i = 0; i < counts[iter->X()][iter->Y()].size(); i++)
        {
            counts[iter->X()][iter->Y()][i] -= expected_photons;
        }
    }

    void PhotonCount::Threshold(const Iterator *iter, int threshold)
    {
        for (int i = 0; i < counts[iter->X()][iter->Y()].size(); i++)
        {
            if (counts[iter->X()][iter->Y()][i] < threshold)
            {
                counts[iter->X()][iter->Y()][i] = 0;
            }
        }
    }

    vector<bool> PhotonCount::AboveThreshold(const Iterator* iter, int threshold) const
    {
        vector<int> data = counts[iter->X()][iter->Y()];
        vector<bool> triggers = vector<bool>();
        for (int i = 0; i < data.size(); i++) triggers.push_back(data[i] > threshold);
        return triggers;
    }

    int PhotonCount::FindThreshold(double global_rate, double sigma)
    {
        double max_prob = gauss.Integral(sigma, Infinity());
        double mean = RealNoiseRate(global_rate);
        int thresh = (int) Floor(sigma * Sqrt(mean));
        while(PoissonSum(mean, thresh) > max_prob) thresh++;
        return thresh;
    }

    void PhotonCount::Subset(vector<vector<vector<bool>>> good_bins)
    {
        for (int i = 0; i < counts.size(); i++)
        {
            for (int j = 0; j < counts[i].size(); j++)
            {
                for (int t = 0; t < counts[i][j].size(); t++)
                {
                    if (!good_bins[i][j][t]) counts[i][j][t] = 0;
                }
            }
        }
    }

    void PhotonCount::Equalize()
    {
        if (equalized) return;
        Iterator iter = GetIterator();
        while (iter.Next()) ExpandVector(iter.X(), iter.Y(), NBins());
        equalized = true;
    }

    bool PhotonCount::ValidPixel(int x_index, int y_index) const
    {
        int n_half = n_pixels / 2;
        double arc_length = n_half * linear_size;
        double angle = n_half * angular_size;
        double radius = arc_length / angle;
        double max_dev = radius * Sin(angle);
        TVector3 direction = Direction(x_index, y_index) * radius;
        return Utility::WithinXYDisk(direction, max_dev);
    }

    TVector3 PhotonCount::Direction(int x_index, int y_index) const
    {
        // Reverse the operations of AddPhoton. Add 0.5 to account for the Floor() operation.
        double pixels_up = (y_index - n_pixels / 2.0 + 0.5);
        double elevation = pixels_up * angular_size;
        double pixels_horiz = (x_index - n_pixels / 2.0 + 0.5);
        double azimuth = pixels_horiz * angular_size * Cos(elevation);

        // A positive azimuth should correspond to a positive x component.
        return TVector3(Cos(elevation) * Sin(azimuth), Sin(elevation), Cos(elevation) * Cos(azimuth));
    }

    void PhotonCount::ExpandVector(int x_index, int y_index, int size)
    {
        if (counts[x_index][y_index].size() < size) counts[x_index][y_index].resize(size, 0);
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
        for (int i = 0; Poisson(mean, i) > 10e-6 * Poisson(mean, min); i++)
        {
            sum += Poisson(mean, i);
        }
        return sum;
    }

    double PhotonCount::Poisson(double mean, int x)
    {
        return Exp(-mean) * Power(x, mean) / Factorial(x);
    }
}

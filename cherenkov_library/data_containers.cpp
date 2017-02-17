// data_containers.cpp
// cherenkov_library
//
// Created by Matthew Dutson on 9/8/16.
//
//

#include "data_containers.h"
#include "TMath.h"
#include "utility.h"

using namespace TMath;
using std::vector;

namespace cherenkov_library
{
    SignalIterator::SignalIterator(std::vector<std::vector<bool>> validPixels)
    {
        this->valid_pixels = validPixels;
        Reset();
    }

    int SignalIterator::X()
    {
        return x_current;
    }

    int SignalIterator::Y()
    {
        return y_current;
    }

    bool SignalIterator::Next()
    {
        // If the outer vector is empty, then there will never be anything to iterate through.
        if (valid_pixels.size() == 0)
        {
            return false;
        }

        // Loop until we find a valid pixel.
        bool valid = false;
        while (!valid)
        {
            // If we're at the end of the current column, try to move to the next column. Otherwise, increment y.
            if (y_current == valid_pixels[x_current].size() - 1)
            {
                // If we're at the end of the last column, return false. Otherwise, move to the next column.
                if (x_current == valid_pixels.size() - 1)
                {
                    return false;
                }
                else
                {
                    x_current++;
                    y_current = -1;
                }
            }
            else
            {
                y_current++;
            }
            if (y_current > -1)
            {
                valid = valid_pixels[x_current][y_current];
            }
        }
        return true;
    }

    void SignalIterator::Reset()
    {
        x_current = 0;
        y_current = -1;
    }

    PhotonCount::PhotonCount(int n_pmt_across, double start_time, double time_bin, double pmt_angular_size,
                             double pmt_linear_size)
    {
        // Copy parameters.
        this->n_pixels = n_pmt_across;
        this->start_time = start_time;
        this->bin_size = time_bin;
        this->pmt_angular_size = pmt_angular_size;
        this->pmt_linear_size = pmt_linear_size;

        // We haven't seen any non-noise photons yet, so the last time is the start time and channels are equalized.
        last_time = start_time;
        equalized = true;

        // Initialize the photon count and valid pixel structures.
        photon_counts = vector<vector<vector<int>>>();
        valid_pixels = vector<vector<bool>>();
        for (int i = 0; i < n_pixels; i++)
        {
            photon_counts.push_back(vector<vector<int>>());
            valid_pixels.push_back(vector<bool>());
            for (int j = 0; j < n_pmt_across; j++)
            {
                photon_counts[i].push_back(vector<int>());
                valid_pixels[i].push_back(ValidPixel(i, j));
            }
        }
    }

    int PhotonCount::Size()
    {
        return n_pixels;
    }

    int PhotonCount::NBins()
    {
        return Bin(last_time) + 1;
    }

    double PhotonCount::Time(int bin)
    {
        return bin * bin_size + start_time;
    }

    int PhotonCount::Bin(double time)
    {
        return Floor((time - start_time) / bin_size);
    }

    void PhotonCount::AddPhoton(double time, TVector3 direction, double thinning)
    {
        // If the time is earlier than the time of the first index, we would run into time binning errors.
        if (time < start_time) return;

        // The angle of the yz projection with the z axis
        double elevation = ATan(direction.Y() / direction.Z());

        // The angle of the xz projection with the z axis
        double azimuth = ATan(direction.X() / direction.Z());

        // The number of photomultipliers across is equal to the size of the data container.
        int y_index = Floor(elevation / pmt_angular_size) + n_pixels / 2;
        int x_index = Floor(azimuth / (pmt_angular_size * Cos(elevation))) + n_pixels / 2;

        // If the location is valid, add the pixel to the underlying data structure.
        if (ValidPixel(x_index, y_index))
        {
            // Make sure the vector is big enough.
            int time_slot = Bin(time);
            ExpandVector(x_index, y_index, time_slot + 1);

            // Add the photon.
            photon_counts[x_index][y_index][time_slot] += thinning;

            // Keep track of the latest photon seen.
            if (time > last_time)
            {
                last_time = time;
            }

            // We may have disrupted the equalization of the channels.
            equalized = false;
        }
    }

    void PhotonCount::EqualizeTimeSeries()
    {
        // If needed, resize all channels so they have bins up through the last photon seen.
        if (equalized)
        {
            return;
        }
        else
        {
            SignalIterator iter = Iterator();
            while (iter.Next())
            {
                ExpandVector(iter.X(), iter.Y(), NBins());
            }
            equalized = true;
        }
    }

    void PhotonCount::AddNoise(double noise_rate, SignalIterator current, TRandom3* rng)
    {
        // Ensure all channels have the same length.
        EqualizeTimeSeries();

        // Randomly determine the number of photons in each bin according to the Poisson distribution.
        double expected_photons = RealNoiseRate(noise_rate);
        for (int i = 0; i < NBins(); i++)
        {
            photon_counts[current.X()][current.Y()][i] += rng->Poisson(expected_photons);
        }
    }

    void PhotonCount::SubtractNoise(double noise_rate, SignalIterator current)
    {
        double expected_photons = RealNoiseRate(noise_rate);
        for (int i = 0; i < NBins(); i++)
        {
            photon_counts[current.X()][current.Y()][i] -= expected_photons;
        }
    }

    double PhotonCount::RealNoiseRate(double universal_rate)
    {
        // Calculate the number of noise photons per second from the number of photons per second per steradian.
        double solid_angle = Sq(pmt_angular_size);
        double pixel_rate = universal_rate * solid_angle;

        // Determine the mean number of Poisson events in a single time bin.
        return bin_size * pixel_rate;
    }

    void PhotonCount::ClearNoise(SignalIterator current, double noise_rate, double hold_thresh)
    {
        // Determine the minimum threshold for non-noise signals.
        double real_noise = RealNoiseRate(noise_rate);
        double thresh = hold_thresh * Sqrt(real_noise);

        // Zero any measurements which are below the noise threshold.
        vector<int> data = photon_counts[current.X()][current.Y()];
        for (int i = 0; i < data.size(); i++)
        {
            if (data[i] < thresh)
            {
                data[i] = 0;
            }
        }
    }

    vector<bool> PhotonCount::FindTriggers(SignalIterator current, double noise_rate, double trigger_thresh)
    {
        // Determine the minimum threshold for triggering.
        double real_noise = RealNoiseRate(noise_rate);
        double thresh = trigger_thresh * Sqrt(real_noise);

        // Find all bins in which the count exceeded the triggering threshold.
        vector<int> data = photon_counts[current.X()][current.Y()];
        vector<bool> triggers = vector<bool>();
        for (int i = 0; i < data.size(); i++)
        {
            if (data[i] > thresh)
            {
                triggers.push_back(true);
            }
            else
            {
                triggers.push_back(false);
            }
        }
        return triggers;
    }

    void PhotonCount::EraseNonTriggered(vector<bool> good_bins)
    {
        SignalIterator iter = Iterator();
        while (iter.Next())
        {
            vector<int> data = photon_counts[iter.X()][iter.Y()];
            for (int i = 0; i < data.size(); i++)
            {
                if (!good_bins.at(i))
                {
                    data[i] = 0;
                }
            }
        }
    }

    TVector3 PhotonCount::Direction(SignalIterator current)
    {
        return Direction(current.X(), current.Y());
    }

    vector<int> PhotonCount::Signal(SignalIterator current)
    {
        return photon_counts[current.X()][current.Y()];
    }

    vector<vector<bool>> PhotonCount::GetValid()
    {
        return valid_pixels;
    }

    SignalIterator PhotonCount::Iterator()
    {
        return SignalIterator(valid_pixels);
    }

    int PhotonCount::SumBins(SignalIterator current)
    {
        vector<int> data = photon_counts[current.X()][current.Y()];
        int sum = 0;
        for (int count : data)
        {
            sum += count;
        }
        return sum;
    }

    double PhotonCount::AverageTime(SignalIterator iter)
    {
        int sum = SumBins(iter);
        double average = 0;
        vector<int> pixel_data = photon_counts[iter.X()][iter.Y()];
        for (int i = 0; i < pixel_data.size(); i++)
        {
            average += Time(i) * pixel_data[i] / (double) sum;
        }
        return average;
    }

    bool PhotonCount::ValidPixel(int x_index, int y_index)
    {
        int n_half = n_pixels / 2;
        double arc_length = n_half * pmt_linear_size;
        double angle = n_half * pmt_angular_size;
        double radius = arc_length / angle;
        double max_dev = radius * Sin(angle);
        TVector3 direction = Direction(x_index, y_index) * radius;
        return WithinXYDisk(direction, max_dev);
    }

    TVector3 PhotonCount::Direction(int x_index, int y_index)
    {
        // Reverse the operations of AddPhoton. Add 0.5 to account for the Floor() operation.
        double pixels_up = (y_index - n_pixels / 2.0 + 0.5);
        double elevation = pixels_up * pmt_angular_size;
        double pixels_horiz = (x_index - n_pixels / 2.0 + 0.5);
        double azimuth = pixels_horiz * pmt_angular_size * Cos(elevation);

        // TODO: Check that these rotations are in the correct direction (perhaps with unit testing).
        // A positive azimuth should correspond to a positive x component.
        return TVector3(Cos(elevation) * Sin(azimuth), Sin(elevation), Cos(elevation) * Cos(azimuth));
    }

    void PhotonCount::ExpandVector(int x_index, int y_index, int min_size)
    {
        if (photon_counts[x_index][y_index].size() < min_size)
        {
            photon_counts[x_index][y_index].resize(min_size, 0);
        }
    }
}

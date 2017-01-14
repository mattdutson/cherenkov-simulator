// data_containers.cpp
// cherenkov_simulator
//
// Created by Matthew Dutson on 9/8/16.
//
//

#include "data_containers.h"
#include "TMath.h"

using namespace TMath;
using std::vector;

namespace cherenkov_simulator
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
        this->start_time = start_time;
        this->bin_size = time_bin;
        this->pmt_angular_size = pmt_angular_size;
        this->pmt_linear_size = pmt_linear_size;

        // We haven't seen any non-noise photons yet, so the last time is the start time.
        last_time = start_time;

        // Initialize the photon count and valid pixel structures.
        photon_counts = vector<vector<vector<int>>>();
        valid_pixels = vector<vector<bool>>();
        for (int i = 0; i < n_pmt_across; i++)
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

    void PhotonCount::AddPhoton(double time, TVector3 direction)
    {
        // If the time is earlier than the time of the first index, we would run into time binning errors.
        if (time < start_time) return;

        // The angle of the yz projection with the z axis
        double elevation = ATan(direction.Y() / direction.Z());

        // The angle of the xz projection with the z axis
        double azimuth = ATan(direction.X() / direction.Z());

        // The number of photomultipliers across is equal to the size of the data container.
        int y_index = Floor(elevation / pmt_angular_size) + photon_counts.size() + 1;
        int x_index = Floor(azimuth / (pmt_angular_size * Cos(elevation))) + photon_counts.size() + 1;

        // If the location is valid, add the pixel to the underlying data structure.
        if (ValidPixel(x_index, y_index))
        {
            // Make sure the vector is big enough.
            int time_slot = Bin(time);
            ExpandVector(x_index, y_index, time_slot + 1);

            // Add the photon.
            photon_counts[x_index][y_index][time_slot]++;

            // Keep track of the latest photon seen.
            if (time > last_time)
            {
                last_time = time;
            }
        }
    }

    void PhotonCount::AddNoise(double noise_rate, SignalIterator current, TRandom3 rng)
    {
        // Resize all channels so they have bins up through the last photon seen.
        SignalIterator iter = Iterator();
        while (iter.Next())
        {
            ExpandVector(iter.X(), iter.Y(), NBins());
        }

        // Calculate the number of noise photons per second from the number of photons per second per steradian per
        // square centimeter.
        double solid_angle = Sq(pmt_angular_size);
        double area = Sq(pmt_linear_size);
        double pixel_rate = noise_rate * solid_angle * area;

        // Determine the mean number of Poisson events in a single time bin.
        double expected_photons = bin_size * pixel_rate;

        // Randomly determine the number of photons in each bin according to the Poisson distribution.
        vector<int> data = photon_counts[current.X()][current.Y()];
        for (int i = 0; i < data.size(); i++)
        {
            data[i] += rng.Poisson(expected_photons);
        }
    }

    TVector3 PhotonCount::Direction(SignalIterator current)
    {
        // Reverse the operations of AddPhoton. Add 0.5 to account for the Floor() operation.
        double elevation = (current.Y() - 1.0 - photon_counts.size() + 0.5) * pmt_angular_size;
        double azimuth = (current.X() - 1.0 - photon_counts.size() + 0.5) * pmt_angular_size * Cos(elevation);

        // TODO: Check that these rotations are in the correct direction (perhaps with unit testing).
        TVector3 direction = TVector3(0, 0, 1);
        direction.RotateX(elevation);
        direction.RotateY(azimuth);
        return direction;
    }

    vector<int> PhotonCount::Signal(SignalIterator current)
    {
        return photon_counts[current.X()][current.Y()];
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

    bool PhotonCount::ValidPixel(int x_index, int y_index)
    {
        int n_pixels = photon_counts.size();
        double center_index = (n_pixels + 1) / 2;
        double x_dist = x_index - center_index;
        double y_dist = y_index - center_index;

        // The radius of the valid circle is equal to half the number of valid pixels.
        return Sq(x_dist) + Sq(y_dist) < Sq(n_pixels / 2.0);
    }

    void PhotonCount::ExpandVector(int x_index, int y_index, int min_size)
    {
        if (photon_counts[x_index][y_index].size() < min_size)
        {
            photon_counts[x_index][y_index].resize(min_size, 0);
        }
    }
}

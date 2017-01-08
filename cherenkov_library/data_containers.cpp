// data_containers.cpp
// cherenkov_simulator
//
// Created by Matthew Dutson on 9/8/16.
//
//

#include "data_containers.h"
#include "TMath.h"

using namespace std;
using namespace TMath;

namespace cherenkov_simulator
{

    PhotonCount::PhotonCount(int n_pmt_across, double start_time, double time_bin, double pmt_angular_size)
    {
        // Copy parameters over.
        this->start_time = start_time;
        this->time_bin = time_bin;
        this->pmt_angular_size = pmt_angular_size;

        // Initialize the data structure.
        photon_counts = vector<vector<vector<int>>>();
        for (int i = 0; i < n_pmt_across; i++)
        {
            photon_counts.push_back(vector<vector<int>>());
            for (int j = 0; j < n_pmt_across; j++)
            {
                photon_counts[i].push_back(vector<int>());
            }
        }
    }

    void PhotonCount::AddPhoton(double time, TVector3 direction)
    {
        direction = -direction;

        // The angle of the yz projection with the z axis. See notes for details.
        double elevation = ATan(direction.Y() / direction.Z());

        // The angle of the xz projection with the z axis. See notes for details.
        double azimuth = ATan(direction.X() / direction.Z());

        // The number of photomultipliers across is equal to the size of the data container.
        int y_index = Floor(elevation / pmt_angular_size) + photon_counts.size() + 1;
        int x_index = Floor(azimuth / (pmt_angular_size * Cos(elevation))) + photon_counts.size() + 1;

        // If the location is valid, add the pixel to the underlying data structure.
        if (!ValidPixel(x_index, y_index))
        {
            return;
        } else
        {
            int time_slot = Floor((time - start_time) / time_bin);
            if (photon_counts[x_index][y_index].size() < time_slot)
            {
                photon_counts[x_index][y_index].resize(time_slot, 0);
            }
            photon_counts[x_index][y_index][time_slot]++;
        }
    }

    bool PhotonCount::ValidPixel(int x_index, int y_index)
    {
        int n_pixels = photon_counts.size();
        double center_index = (n_pixels + 1) / 2;
        double x_dist = x_index - center_index;
        double y_dist = y_index - center_index;
        if (x_dist * x_dist + y_dist * y_dist > n_pixels * n_pixels / 4.0)
        {
            return false;
        } else
        {
            return true;
        }
    }
}

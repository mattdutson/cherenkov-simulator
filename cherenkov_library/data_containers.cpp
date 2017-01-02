// data_containers.cpp
// cherenkov_simulator
//
// Created by Matthew Dutson on 9/8/16.
//
//

#include "data_containers.h"

using namespace std;

namespace cherenkov_simulator
{

    PhotonCount::PhotonCount(int pixel_count)
    {
        photon_counts = vector<vector<vector<int>>>();
        for (int i = 0; i < pixel_count; i++)
        {
            photon_counts.push_back(vector<vector<int>>());
            for (int j = 0; j < pixel_count; j++)
            {
                photon_counts[i].push_back(vector<int>());
            }
        }
    }

    void PhotonCount::AddPhoton(int time_slot, int x_index, int y_index)
    {
        if (!ValidPixel(x_index, y_index))
        {
            return;
        } else
        {
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

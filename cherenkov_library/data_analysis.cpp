//
// Created by Matthew Dutson on 1/9/17.
//

#include "data_analysis.h"
#include <array>

using std::array;
using std::vector;

namespace cherenkov_simulator
{
    bool Above(TVector3 reference, TVector3 other)
    {
        if (other.Z() > reference.Z())
        {
            return true;
        }
        else if (other.Z() < reference.Z())
        {
            return false;
        }
        else
        {
            if (other.Y() > reference.Y())
            {
                return true;
            }
            else
            {
                return false;
            }
        }
    }

    vector<array<double, 2>> CollapseToProfile(PhotonCount data, Plane s_d_plane, TVector3 shower_axis)
    {
        // TODO: Adjust photon count data container so each channel has a vector of the same length.
        vector<array<double, 2>> collapsed_bins = vector<array<double, 2>>();
        SignalIterator iter = data.Iterator();
        TVector3 norm = s_d_plane.Normal();
        TVector3 axis_project = shower_axis - shower_axis.Dot(norm) * norm;
        while (iter.Next())
        {
            TVector3 direction = data.Direction(iter);
            TVector3 projection = direction - direction.Dot(norm) * norm;
            double angle;
            if (Above(axis_project, projection))
            {
                angle = projection.Angle(axis_project);
            }
            else
            {
                angle = -projection.Angle(axis_project);
            }
            array<double, 2> point = {angle, (double) data.SumBins(iter)};
            collapsed_bins.push_back(point);
        }

        // TODO: Figure out a way to bin photon counts in a way that doesn't lead to aliasing.
        return collapsed_bins;
    }

    vector<array<double, 2>> SuperimposeTimes(PhotonCount data)
    {
        // Initialize the structure and find x-axis labels (times).
        vector<array<double, 2>> superposition = vector<array<double, 2>>();
        for (int i = 0; i < data.NBins(); i++)
        {
            superposition.push_back({data.Time(i), 0.0});
        }

        // Iterate over all time signals.
        SignalIterator iter = data.Iterator();
        while (iter.Next())
        {
            vector<int> signal = data.Signal(iter);
            for (int i = 0; i < signal.size(); i++)
            {
                superposition[i][1] += signal[i];
            }
        }
        return superposition;
    };
}
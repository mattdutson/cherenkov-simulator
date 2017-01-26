//
// Created by Matthew Dutson on 1/9/17.
//

#include "data_analysis.h"
#include <array>

using std::array;
using std::vector;

namespace cherenkov_library
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

    vector<vector<double>> CollapseToProfile(PhotonCount data, Plane s_d_plane, TVector3 shower_axis)
    {
        // TODO: Adjust photon count data container so each channel has a vector of the same length.
        vector<vector<double>> collapsed_bins = vector<vector<double>>();
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
            vector<double> point = {angle, (double) data.SumBins(iter)};
            collapsed_bins.push_back(point);
        }

        // TODO: Figure out a way to bin photon counts in a way that doesn't lead to aliasing.
        return collapsed_bins;
    }

    vector<vector<double>> SuperimposeTimes(PhotonCount data)
    {
        // Initialize the structure and find x-axis labels (times).
        vector<vector<double>> superposition = vector<vector<double>>();
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

    TH2C GetValidMap(PhotonCount data)
    {
        std::vector<std::vector<bool>> valid = data.GetValid();
        TH2C histo = TH2C("Valid", "Valid", valid.size(), 0, valid.size(), valid[0].size(), 0, valid[0].size());
        for (int i = 0; i < valid.size(); i++)
        {
            for (int j = 0; j < valid[i].size(); j++)
            {
                if (valid[i][j])
                {
                    // The zeroth bin is the underflow, so start at 1.
                    histo.SetBinContent(i + 1, j + 1, 1.0);
                }
                else
                {
                    // The zeroth bin is the underflow, so start at 1.
                    histo.SetBinContent(i + 1, j + 1, 0.0);
                }
            }
        }
        return histo;
    }
}
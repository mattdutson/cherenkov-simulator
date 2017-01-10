//
// Created by Matthew Dutson on 1/9/17.
//

#include "analysis.h"

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

    std::vector<int[2]> CollapseToProfile(PhotonCount data, Plane s_d_plane, TVector3 shower_axis)
    {
        // TODO: Adjust photon count data container so each channel has a vector of the same length.
        std::vector<int[2]> collapsed_bins = std::vector<int[2]>();
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
            collapsed_bins.push_back({angle, data.SumBins(iter)});
        }

        // TODO: Figure out a way to bin photon counts in a way that doesn't lead to aliasing.
        return collapsed_bins;
    }
}
// Analysis.cpp
//
// Author: Matthew Dutson
//
// Implementation of Analysis.h

#include "Analysis.h"
#include "Utility.h"

using std::vector;

namespace cherenkov_lib
{
    void
    Analysis::CollapseToProfile(PhotonCount data, Plane s_d_plane, TVector3 shower_axis, vector<double>* angles,
                                vector<double>* counts)
    {
        // TODO: Figure out a way to bin photon counts in a way that doesn't lead to aliasing.
        *angles = vector<double>();
        *counts = vector<double>();
        SignalIterator iter = data.Iterator();
        TVector3 norm = s_d_plane.Normal();
        TVector3 axis_project = shower_axis - shower_axis.Dot(norm) * norm;
        while (iter.Next())
        {
            TVector3 direction = data.Direction(&iter);
            TVector3 projection = direction - direction.Dot(norm) * norm;
            double angle;
            if (Utility::Above(axis_project, projection))
            {
                angle = projection.Angle(axis_project);
            }
            else
            {
                angle = -projection.Angle(axis_project);
            }
            angles->push_back(angle);
            counts->push_back(data.SumBins(&iter));
        }
    }

    void Analysis::SuperimposeTimes(PhotonCount data, vector<double>* times, vector<double>* counts)
    {
        // Initialize the structure and find x-axis labels (times).
        *times = vector<double>();
        *counts = vector<double>();
        for (int i = 0; i < data.NBins(); i++)
        {
            times->push_back(data.Time(i));
            counts->push_back(0.0);
        }

        // Iterate over all pixels and add their time signals.
        SignalIterator iter = data.Iterator();
        while (iter.Next())
        {
            vector<int> signal = data.Signal(&iter);
            for (int i = 0; i < signal.size(); i++)
            {
                counts->at(i) += signal[i];
            }
        }
    };

    TH2C Analysis::GetValidMap(PhotonCount data)
    {
        return GetBooleanMap(data.GetValid());
    }

    TH2C Analysis::GetBooleanMap(vector<vector<bool>> valid)
    {
        TH2C histo = TH2C("Valid", "Valid", valid.size(), 0, valid.size(), valid[0].size(), 0, valid[0].size());
        for (int i = 0; i < valid.size(); i++)
        {
            for (int j = 0; j < valid[i].size(); j++)
            {
                // The zeroth bin is the underflow, so start at 1.
                histo.SetBinContent(i + 1, j + 1, valid[i][j] ? 1.0 : 0.0);

            }
        }
        return histo;
    }

    TGraph Analysis::MakeProfileGraph(PhotonCount data)
    {
        vector<double> times, counts;
        SuperimposeTimes(data, &times, &counts);
        return TGraph(times.size(), &(times[0]), &(counts[0]));
    }

    TH2I Analysis::MakeSumMap(PhotonCount data)
    {
        int size = data.Size();
        TH2I histo = TH2I("Sums", "Sums", size, 0, size, size, 0, size);
        SignalIterator iter = data.Iterator();
        while (iter.Next())
        {
            // The underflow seems to be defined differently for the TH1I than for the TH1C
            histo.Fill(iter.X(), iter.Y(), 1.0);
        }
        return histo;
    }
}
// Analysis.cpp
//
// Author: Matthew Dutson
//
// Implementation of Analysis.h

#include "Analysis.h"

using namespace std;

namespace cherenkov_simulator
{
    TGraph Analysis::MakeTimeProfile(const PhotonCount& data)
    {
        Double1D times, counts;
        SuperimposeTimes(data, times, counts);
        TGraph output = TGraph((int) times.size(), &(times[0]), &(counts[0]));
        output.SetTitle("Detector Time Profile");
        output.GetXaxis()->SetTitle("Time (s)");
        output.GetYaxis()->SetTitle("Total Photons Seen");
        return output;
    }

    TH2I Analysis::MakePixlProfile(const PhotonCount& data, string name, bool reverse_y)
    {
        auto size = (int) data.Size();
        TH2I histo = TH2I(name.c_str(), "Bin Signal Sums", size, 0, size, size, 0, size);
        histo.SetXTitle("x Bin");
        histo.SetYTitle("y Bin");
        PhotonCount::Iterator iter = data.GetIterator();
        while (iter.Next())
        {
            int y = iter.Y();
            if (reverse_y) y = size - y - 1;
            histo.Fill(iter.X(), y, data.SumBins(iter));
        }
        return histo;
    }

    TH2C Analysis::GetValidMap(const PhotonCount& data)
    {
        return GetBooleanMap(data.GetValid());
    }

    TH2C Analysis::GetBooleanMap(const Bool2D& valid)
    {
        TH2C histo = TH2C("valid_map", "Valid Pixels", (int) valid.size(), 0, valid.size(), (int) valid[0].size(), 0, valid[0].size());
        histo.SetXTitle("x Bin");
        histo.SetYTitle("y Bin");

        // The zeroth bin is the underflow, so start at 1.
        for (int i = 0; i < valid.size(); i++)
            for (int j = 0; j < valid[i].size(); j++)
                histo.SetBinContent(i + 1, j + 1, valid[i][j] ? 1.0 : 0.0);
        return histo;
    }

    void Analysis::SuperimposeTimes(const PhotonCount& data, Double1D& times, Double1D& counts)
    {
        times = Double1D();
        counts = Double1D();
        for (int i = 0; i < data.NBins(); i++)
        {
            times.push_back(data.Time(i));
            counts.push_back(0.0);
        }

        PhotonCount::Iterator iter = data.GetIterator();
        while (iter.Next())
        {
            Short1D signal = data.Signal(iter);
            for (size_t i = 0; i < signal.size(); i++)
                counts.at(i) += signal[i];
        }
    }
}
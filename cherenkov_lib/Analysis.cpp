// Analysis.cpp
//
// Author: Matthew Dutson
//
// Implementation of Analysis.h

#include "Analysis.h"

using std::vector;

namespace cherenkov_simulator
{
    TGraph Analysis::MakeProfileGraph(const PhotonCount& data)
    {
        Double1D times, counts;
        SuperimposeTimes(data, times, counts);
        TGraph output = TGraph((int) times.size(), &(times[0]), &(counts[0]));
        output.SetTitle("Detector Time Profile");
        output.GetXaxis()->SetTitle("Time (s)");
        output.GetYaxis()->SetTitle("Total Photons Seen");
        return output;
    }

    TH2I Analysis::MakeSumMap(const PhotonCount& data, bool reverse_y)
    {
        int size = (int) data.Size();
        TH2I histo = TH2I("sum_map", "Bin Signal Sums", size, 0, size, size, 0, size);
        histo.SetXTitle("x Bin");
        histo.SetYTitle("y Bin");
        PhotonCount::Iterator iter = data.GetIterator();
        while (iter.Next())
        {
            size_t y = iter.Y();
            if (reverse_y) y = size - 1 - y;
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
        TH2C histo = TH2C("valid_map", "Valid Pixels", (int) valid.size(), 0, valid.size(), (int) valid[0].size(), 0,
                          valid[0].size());
        histo.SetXTitle("x Bin");
        histo.SetYTitle("y Bin");
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

    void Analysis::SuperimposeTimes(const PhotonCount& data, Double1D& times, Double1D& counts)
    {
        // Initialize the structure and find x-axis labels (times).
        times = Double1D();
        counts = Double1D();
        for (int i = 0; i < data.NBins(); i++)
        {
            times.push_back(data.Time(i));
            counts.push_back(0.0);
        }

        // Iterate over all pixels and add their time signals.
        PhotonCount::Iterator iter = data.GetIterator();
        while (iter.Next())
        {
            Int1D signal = data.Signal(iter);
            for (size_t i = 0; i < signal.size(); i++)
            {
                counts.at(i) += signal[i];
            }
        }
    };
}
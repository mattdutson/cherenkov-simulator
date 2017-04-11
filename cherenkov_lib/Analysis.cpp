// Analysis.cpp
//
// Author: Matthew Dutson
//
// Implementation of Analysis.h

#include "Analysis.h"
#include "Utility.h"

using std::vector;

namespace cherenkov_simulator
{
    TGraph Analysis::MakeProfileGraph(const PhotonCount& data)
    {
        vector<double> times, counts;
        SuperimposeTimes(data, times, counts);
        TGraph output = TGraph(times.size(), &(times[0]), &(counts[0]));
        output.SetTitle("Detector Time Profile");
        output.GetXaxis()->SetTitle("Time (s)");
        output.GetYaxis()->SetTitle("Total Photons Seen");
        return output;
    }

    TH2I Analysis::MakeSumMap(const PhotonCount& data)
    {
        int size = data.Size();
        TH2I histo = TH2I("sum_map", "Bin Signal Sums", size, 0, size, size, 0, size);
        histo.SetXTitle("x Bin");
        histo.SetYTitle("y Bin");
        PhotonCount::Iterator iter = data.GetIterator();
        while (iter.Next())
        {
            // The underflow seems to be defined differently for the TH1I than for the TH1C
            histo.Fill(iter.X(), iter.Y(), data.SumBins(iter));
        }
        return histo;
    }

    TH2C Analysis::GetValidMap(const PhotonCount& data)
    {
        return GetBooleanMap(data.GetValid());
    }

    TH2C Analysis::GetBooleanMap(const vector<vector<bool>>& valid)
    {
        TH2C histo = TH2C("valid_map", "Valid Pixels", valid.size(), 0, valid.size(), valid[0].size(), 0, valid[0].size());
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

    void Analysis::SuperimposeTimes(const PhotonCount& data, vector<double>& times, vector<double>& counts)
    {
        // Initialize the structure and find x-axis labels (times).
        times = vector<double>();
        counts = vector<double>();
        for (int i = 0; i < data.NBins(); i++)
        {
            times.push_back(data.Time(i));
            counts.push_back(0.0);
        }

        // Iterate over all pixels and add their time signals.
        PhotonCount::Iterator iter = data.GetIterator();
        while (iter.Next())
        {
            vector<int> signal = data.Signal(iter);
            for (int i = 0; i < signal.size(); i++)
            {
                counts.at(i) += signal[i];
            }
        }
    };
}
// DataStructures.cpp
//
// Author: Matthew Dutson
//
// Implementation of DataStructures.h.

#include <TMath.h>
#include <TRandom3.h>

#include "DataStructures.h"

using namespace std;
using namespace TMath;

namespace cherenkov_simulator
{
    PhotonCount::Iterator::Iterator(Bool2D validPixels)
    {
        this->valid = move(validPixels);
        Reset();
    }

    int PhotonCount::Iterator::X() const
    {
        if (curr_y < 0)
            throw out_of_range("Call Next() before checking the iterator position");
        return curr_x;
    }

    int PhotonCount::Iterator::Y() const
    {
        if (curr_y < 0)
            throw out_of_range("Call Next() before checking the iterator position");
        return curr_y;
    }

    bool PhotonCount::Iterator::Next()
    {
        if (valid.empty()) return false;

        bool found = false;
        while (!found)
        {
            // If we're at the end of the current column, try to move to the next.
            if (curr_y == valid[curr_x].size() - 1)
            {
                if (curr_x == valid.size() - 1) return false;
                curr_x++;
                curr_y = -1;
            }
            else
            {
                curr_y++;
            }
            if (curr_y > -1) found = valid[curr_x][curr_y];
        }
        return true;
    }

    void PhotonCount::Iterator::Reset()
    {
        curr_x = 0;
        curr_y = -1;
    }

    PhotonCount::PhotonCount()
    {
        n_pixels = 0;
        min_time = 0;
        max_time = 0;
        empty = true;
    }

    PhotonCount::PhotonCount(Params params, double min_time, double max_time)
    {
        n_pixels = params.n_pixels;
        bin_size = params.bin_size;
        ang_size = params.ang_size;
        lin_size = params.lin_size;
        this->min_time = min_time;
        this->max_time = max_time;

        trimd = false;
        empty = true;
        frst_time = max_time;
        last_time = min_time;

        if (n_pixels % 2 != 0)
            throw invalid_argument("Number of pixels must be even");
        if (bin_size <= 0.0)
            throw invalid_argument("Bin size must be positive");
        if (ang_size <= 0.0)
            throw invalid_argument("Angular size must be positive");
        if (lin_size <= 0.0)
            throw invalid_argument("Linear size must be positive");
        if (NBins() * Sq(n_pixels) * sizeof(short) > params.max_byte)
            throw out_of_range("Warning: too much memory requested due to shower direction");

        counts = Short3D(n_pixels, Short2D(n_pixels, Short1D(NBins(), 0)));
        sums = Short2D(n_pixels, Short1D(n_pixels, 0));
        valid = Bool2D(n_pixels, Bool1D(n_pixels, false));
        for (int i = 0; i < NBins(); i++)
            for (int j = 0; j < NBins(); j++)
                valid[i][j] = IsValid(i, j);
    }

    Bool2D PhotonCount::GetValid() const
    {
        return valid;
    }

    size_t PhotonCount::Size() const
    {
        return n_pixels;
    }

    size_t PhotonCount::NBins() const
    {
        return min_time == max_time ? 0 : Bin(max_time) + 1;
    }

    bool PhotonCount::Empty() const
    {
        return empty;
    }

    double PhotonCount::Time(int bin) const
    {
        if (bin < 0 || bin >= NBins())
            throw out_of_range("Invalid bin");
        return bin * bin_size + min_time + bin_size / 2.0;
    }

    size_t PhotonCount::Bin(double time) const
    {
        if (time < min_time || time > max_time)
            throw out_of_range("Invalid time");
        return (size_t) Floor((time - min_time) / bin_size);
    }

    double PhotonCount::DetectorAxisAngle() const
    {
        return n_pixels / 2.0 * ang_size;
    }

    TVector3 PhotonCount::Direction(const Iterator& iter) const
    {
        return Direction(iter.X(), iter.Y());
    }

    Short1D PhotonCount::Signal(const Iterator& iter) const
    {
        return counts[iter.X()][iter.Y()];
    }

    int PhotonCount::SumBins(const Iterator& iter) const
    {
        return sums[iter.X()][iter.Y()];
    }

    int PhotonCount::SumBinsFiltered(const Iterator& iter, const Bool3D& filter) const
    {
        int sum = 0;
        for (size_t i = 0; i < NBins(); i++)
            sum += filter[iter.X()][iter.Y()][i] ? counts[iter.X()][iter.Y()][i] : 0;
        return sum;
    }

    double PhotonCount::AverageTime(const Iterator& iter) const
    {
        int sum = SumBins(iter);
        if (sum == 0)
            throw invalid_argument("Channel is empty, division by zero");
        double average = 0;
        for (int i = 0; i < NBins(); i++)
            average += counts[iter.X()][iter.Y()][i] * Time(i) / sum;
        return average;
    }

    double PhotonCount::TimeError(const Iterator& iter) const
    {
        int sum = SumBins(iter);
        double mean = AverageTime(iter);
        double variance = 0;
        for (int i = 0; i < NBins(); i++)
            variance += counts[iter.X()][iter.Y()][i] * Sq(Time(i) - mean) / sum;

        // Add a Sheppard correction before computing the standard deviation.
        variance += Sq(bin_size) / 12.0;
        return Sqrt(variance / sum);
    }

    PhotonCount::Iterator PhotonCount::GetIterator() const
    {
        return Iterator(valid);
    }

    Bool3D PhotonCount::GetFalseMatrix() const
    {
        return Bool3D(Size(), Bool2D(Size(), Bool1D(NBins(), false)));
    }

    void PhotonCount::AddPhoton(double time, TVector3 position, int thinning)
    {
        if (position == TVector3())
            throw invalid_argument("Direction cannot be a zero vector");
        if (time < min_time || time > max_time) return;

        TVector3 direction = -position;
        double elevate = ATan2(direction.Y(), direction.Z());
        double azimuth = ATan2(direction.X(), direction.Z());
        auto y_index = (int) (Floor(elevate / ang_size) + n_pixels / 2);
        auto x_index = (int) (Floor(azimuth / ang_size / Cos(elevate)) + n_pixels / 2);
        
        if (IsValid(x_index, y_index))
        {
            IncrementCell(thinning, (size_t) x_index, (size_t) y_index, Bin(time));
            if (time > last_time) last_time = time;
            if (time < frst_time) frst_time = time;
            trimd = false;
        }
    }

    void PhotonCount::AddNoise(double noise_rate, const Iterator& iter)
    {
        Trim();
        double mean = RealNoiseRate(noise_rate);
        for (size_t i = 0; i < NBins(); i++)
            IncrementCell(gRandom->Poisson(mean), iter, i);
    }

    void PhotonCount::Subtract(double noise_rate, const Iterator& iter)
    {
        auto mean = (int) RealNoiseRate(noise_rate);
        for (size_t i = 0; i < NBins(); i++)
            IncrementCell(-mean, iter, i);
    }

    Bool1D PhotonCount::AboveThreshold(const Iterator& iter, int threshold) const
    {
        Bool1D above = Bool1D(NBins());
        for (size_t i = 0; i < NBins(); i++)
            above[i] = counts[iter.X()][iter.Y()][i] > threshold;
        return above;
    }

    int PhotonCount::FindThreshold(double noise_rate, double sigma) const
    {
        double max_prob = Erfc(sigma / Sqrt(2)) / 2.0;
        double mean = RealNoiseRate(noise_rate);
        auto thresh = (int) Floor(sigma * Sqrt(mean));
        while (PoissonSum(mean, thresh) > max_prob)
            thresh++;
        return thresh - 1;
    }

    void PhotonCount::Subset(const Bool3D& good_bins)
    {
        for (size_t i = 0; i < Size(); i++)
            for (size_t j = 0; j < Size(); j++)
                for (size_t t = 0; t < NBins(); t++)
                    if (!good_bins[i][j][t])
                        IncrementCell(-counts[i][j][t], i, j, t);
    }

    void PhotonCount::Trim()
    {
        if (trimd || empty) return;
        for (int i = 0; i < Size(); i++)
        {
            for (int j = 0; j < Size(); j++)
            {
                auto bgn = counts[i][j].begin();
                auto end = counts[i][j].end();
                counts[i][j].erase(bgn + Bin(last_time) + 1, end);
                counts[i][j].erase(bgn, bgn + Bin(frst_time));
            }
        }
        min_time = min_time + Floor((frst_time - min_time) / bin_size) * bin_size;
        max_time = last_time;
        trimd = true;
    }

    void PhotonCount::IncrementCell(int inc, const Iterator& iter, size_t t)
    {
        IncrementCell(inc, (size_t) iter.X(), (size_t) iter.Y(), t);
    }

    void PhotonCount::IncrementCell(int inc, size_t x_index, size_t y_index, size_t t)
    {
        if (inc > 0) empty = false;
        counts[x_index][y_index][t] += inc;
        sums[x_index][y_index] += inc;
    }

    bool PhotonCount::IsValid(int x_index, int y_index) const
    {
        bool in_range = x_index >= 0 && y_index >= 0 && x_index < n_pixels && y_index < n_pixels;
        if (!in_range) return false;
        
        double arc = n_pixels * lin_size / 2.0;
        double ang = n_pixels * ang_size / 2.0;
        double rad = arc / ang;
        TVector3 direction = Direction(x_index, y_index) * rad;
        return Utility::WithinXYDisk(direction, rad * Sin(ang));
    }

    TVector3 PhotonCount::Direction(int x_index, int y_index) const
    {
        double pixels_vert = (y_index - n_pixels / 2.0 + 0.5);
        double pixels_horz = (x_index - n_pixels / 2.0 + 0.5);
        double elevate = pixels_vert * ang_size;
        double azimuth = pixels_horz * ang_size * Cos(elevate);

        // A positive azimuth should correspond to a positive x component.
        return TVector3(Cos(elevate) * Sin(azimuth), Sin(elevate), Cos(elevate) * Cos(azimuth));
    }

    double PhotonCount::RealNoiseRate(double noise_rate) const
    {
        double pixel_rate = noise_rate * Sq(ang_size);
        return bin_size * pixel_rate;
    }

    double PhotonCount::PoissonSum(double mean, int min)
    {
        double sum = 1.0;
        for (int i = 0; i < min; i++)
            sum -= Poisson(mean, i);
        return sum;
    }

    double PhotonCount::Poisson(double mean, int x)
    {
        return Exp(-mean) * Power(mean, x) / Factorial(x);
    }
}

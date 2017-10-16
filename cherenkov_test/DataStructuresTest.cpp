// DataStructuresTest.cpp
//
// Author: Matthew Dutson
//
// Tests of DataStructures.h

#include <gtest/gtest.h>
#include <TFile.h>

#include "DataStructures.h"
#include "Analysis.h"
#include "Helper.h"

using namespace std;
using namespace TMath;

namespace cherenkov_simulator
{
    /*
     * Note: this class will be able to access private member of the PhotonCount and Iterator classes
     */
    class DataStructuresTest : public ::testing::Test
    {
    private:

        PhotonCount::Params test_params;
        PhotonCount empty_data;
        PhotonCount sample_data;

        virtual void SetUp()
        {
            test_params = PhotonCount::Params();
            test_params.n_pixels = 4;
            test_params.max_byte = 4000000;
            test_params.bin_size = 0.1;
            test_params.ang_size = 0.08;
            test_params.lin_size = 2.5;

            empty_data = PhotonCount(test_params, 0.0, 0.95);

            sample_data = PhotonCount(test_params, 0.0, 0.95);
            sample_data.IncrementCell(1, 0, 2, 7);
            sample_data.IncrementCell(2, 1, 1, 3);
            sample_data.IncrementCell(3, 1, 1, 3);
            sample_data.IncrementCell(1, 1, 1, 4);
            sample_data.IncrementCell(8, 1, 1, 9);
            sample_data.trimd = false;
            sample_data.frst_time = 0.35;
            sample_data.last_time = 0.95;
        }

    public:

        PhotonCount CopyEmpty()
        {
            return empty_data;
        }

        PhotonCount CopySample()
        {
            return sample_data;
        }

        PhotonCount::Params CopyParams()
        {
            return test_params;
        }

        double FriendRealNoiseRate(PhotonCount& data, double rate)
        {
            return data.RealNoiseRate(rate);
        }
    };

    /*
     * Test the default constructor.
     */
    TEST_F(DataStructuresTest, DefaultConstruct)
    {
        PhotonCount data = PhotonCount();
        ASSERT_EQ(data.Size(), 0);
        ASSERT_EQ(data.NBins(), 0);
        ASSERT_TRUE(data.Empty());
    }

    /*
     * Test the non-default constructor.
     */
    TEST_F(DataStructuresTest, UserConstructor)
    {
        PhotonCount data = CopyEmpty();
        ASSERT_EQ(4, data.Size());
        ASSERT_EQ(10, data.NBins());
        ASSERT_TRUE(data.Empty());
    }

    /*
     * An invalid_argument exception should be thrown if the number of bins is odd.
     */
    TEST_F(DataStructuresTest, OddNPixels)
    {
        PhotonCount::Params params = CopyParams();
        params.n_pixels = 5;
        try
        {
            PhotonCount(params, 0.0, 0.95);
            FAIL() << "Exception not thrown";
        }
        catch(invalid_argument& err)
        {
            ASSERT_EQ(string("Number of pixels must be even"), err.what());
        }
    }

    /*
     * An invalid_argument exception should be thrown if the size of the underlying data structure is too large.
     */
    TEST_F(DataStructuresTest, OverMaxBytes)
    {
        PhotonCount::Params params = CopyParams();
        try
        {
            params.max_byte = 10;
            PhotonCount(params, 0.0, 0.95);
            FAIL() << "Exception not thrown";
        }
        catch(out_of_range& err)
        {
            ASSERT_EQ(string("Warning: too much memory requested due to shower direction"), err.what());
        }
        params.max_byte = 320;
        PhotonCount(params, 0.0, 0.95);
    }

    /*
     * See what happens if the number of pixels passed to the constructor is zero. This shouldn't cause a division by
     * zero or weird behavior, so it should pass through.
     */
    TEST_F(DataStructuresTest, ZeroPixels)
    {
        PhotonCount::Params params = CopyParams();
        params.n_pixels = 0;
        PhotonCount(params, 0.0, 0.95);
    }

    /*
     * An invalid_argument exception should be thrown if the bin size is non-positive.
     */
    TEST_F(DataStructuresTest, NonPositiveBinSize)
    {
        PhotonCount::Params params = CopyParams();
        params.bin_size = 0.0;
        try
        {
            PhotonCount(params, 0.0, 0.95);
            FAIL() << "Exception not thrown";
        }
        catch(invalid_argument& err)
        {
            ASSERT_EQ(string("Bin size must be positive"), err.what());
        }
        params.bin_size = -0.1;
        try
        {
            PhotonCount(params, 0.0, 0.95);
            FAIL() << "Exception not thrown";
        }
        catch(invalid_argument& err)
        {
            ASSERT_EQ(string("Bin size must be positive"), err.what());
        }
    }

    /*
     * An invalid_argument exception should be thrown if the angular size is non-positive.
     */
    TEST_F(DataStructuresTest, NonPositiveAngularSize)
    {
        PhotonCount::Params params = CopyParams();
        params.ang_size = 0.0;
        try
        {
            PhotonCount(params, 0.0, 0.95);
            FAIL() << "Exception not thrown";
        }
        catch(invalid_argument& err)
        {
            ASSERT_EQ(string("Angular size must be positive"), err.what());
        }
        params.ang_size = -0.08;
        try
        {
            PhotonCount(params, 0.0, 0.95);
            FAIL() << "Exception not thrown";
        }
        catch(invalid_argument& err)
        {
            ASSERT_EQ(string("Angular size must be positive"), err.what());
        }
    }

    /*
     * An invalid_argument exception should be thrown if the linear size is non-positive.
     */
    TEST_F(DataStructuresTest, NonPositiveLinearSize)
    {
        PhotonCount::Params params = CopyParams();
        params.lin_size = 0.0;
        try
        {
            PhotonCount(params, 0.0, 0.95);
            FAIL() << "Exception not thrown";
        }
        catch(invalid_argument& err)
        {
            ASSERT_EQ(string("Linear size must be positive"), err.what());
        }
        params.lin_size = -0.08;
        try
        {
            PhotonCount(params, 0.0, 0.95);
            FAIL() << "Exception not thrown";
        }
        catch(invalid_argument& err)
        {
            ASSERT_EQ(string("Linear size must be positive"), err.what());
        }
    }

    /*
     * Check that the correct pixels are marked as positive.
     */
    TEST_F(DataStructuresTest, ValidPixels)
    {
        PhotonCount data = CopyEmpty();
        Bool2D valid = data.GetValid();
        for (int i = 0; i < valid.size(); i++)
        {
            for (int j = 0; j < valid.size(); j++)
            {
                if (Abs(i - 1.5) < 1 || Abs(j - 1.5) < 1)
                    ASSERT_TRUE(valid[i][j]);
                else
                    ASSERT_FALSE(valid[i][j]);
            }
        }
    }

    /*
     * An out_of_range exception should be thrown if Next() has not been called on a PhotonCount::Iterator before X() or
     * Y() are called.
     */
    TEST_F(DataStructuresTest, InvalidIteratorPosition)
    {
        PhotonCount::Iterator iter = CopyEmpty().GetIterator();
        try
        {
            iter.X();
            FAIL() << "Exception not thrown";
        }
        catch(out_of_range& err)
        {
            ASSERT_EQ(string("Call Next() before checking the iterator position"), err.what());
        }
        try
        {
            iter.Y();
            FAIL() << "Exception not thrown";
        }
        catch(out_of_range& err)
        {
            ASSERT_EQ(string("Call Next() before checking the iterator position"), err.what());
        }
    }

    /*
     * Test the PhotonCount::Iterator. Steps through y and then steps through x. An exception should be thrown if the
     * iterator is at an invalid position and X() or Y() are called.
     */
    TEST_F(DataStructuresTest, TestIterator)
    {
        PhotonCount::Iterator iter = CopyEmpty().GetIterator();

        ASSERT_TRUE(iter.Next());
        ASSERT_EQ(0, iter.X());
        ASSERT_EQ(1, iter.Y());

        ASSERT_TRUE(iter.Next());
        ASSERT_EQ(0, iter.X());
        ASSERT_EQ(2, iter.Y());

        ASSERT_TRUE(iter.Next());
        ASSERT_EQ(1, iter.X());
        ASSERT_EQ(0, iter.Y());

        ASSERT_TRUE(iter.Next());
        ASSERT_EQ(1, iter.X());
        ASSERT_EQ(1, iter.Y());

        ASSERT_TRUE(iter.Next());
        ASSERT_EQ(1, iter.X());
        ASSERT_EQ(2, iter.Y());

        ASSERT_TRUE(iter.Next());
        ASSERT_EQ(1, iter.X());
        ASSERT_EQ(3, iter.Y());

        ASSERT_TRUE(iter.Next());
        ASSERT_EQ(2, iter.X());
        ASSERT_EQ(0, iter.Y());

        ASSERT_TRUE(iter.Next());
        ASSERT_EQ(2, iter.X());
        ASSERT_EQ(1, iter.Y());

        ASSERT_TRUE(iter.Next());
        ASSERT_EQ(2, iter.X());
        ASSERT_EQ(2, iter.Y());

        ASSERT_TRUE(iter.Next());
        ASSERT_EQ(2, iter.X());
        ASSERT_EQ(3, iter.Y());

        ASSERT_TRUE(iter.Next());
        ASSERT_EQ(3, iter.X());
        ASSERT_EQ(1, iter.Y());

        ASSERT_TRUE(iter.Next());
        ASSERT_EQ(3, iter.X());
        ASSERT_EQ(2, iter.Y());

        ASSERT_FALSE(iter.Next());
    }

    /*
     * Ensure that the PhotonCount iterator resets correctly.
     */
    TEST_F(DataStructuresTest, TestIteratorReset)
    {
        PhotonCount::Iterator iter = CopyEmpty().GetIterator();

        iter.Next();
        iter.Next();
        iter.Reset();
        try
        {
            iter.X();
            FAIL() << "Exception not thrown";
        }
        catch(out_of_range& err)
        {
            ASSERT_EQ(string("Call Next() before checking the iterator position"), err.what());
        }
        try
        {
            iter.Y();
            FAIL() << "Exception not thrown";
        }
        catch(exception& err)
        {
            ASSERT_EQ(string("Call Next() before checking the iterator position"), err.what());
        }

        iter.Next();
        ASSERT_EQ(0, iter.X());
        ASSERT_EQ(1, iter.Y());

        iter.Next();
        ASSERT_EQ(0, iter.X());
        ASSERT_EQ(2, iter.Y());
    }

    /*
     * Make sure the empty flag is updated correctly. The container should remain empty if the photon is outside valid
     * time or space bounds. Calling the AddNoise function causes emptiness to be false.
     */
    TEST_F(DataStructuresTest, TestEmpty)
    {
        PhotonCount data0 = CopyEmpty();
        ASSERT_TRUE(data0.Empty());
        data0.AddPhoton(0.2, TVector3(0.0, 0.0, -1.0), 1);
        ASSERT_FALSE(data0.Empty());

        PhotonCount data1 = CopyEmpty();
        ASSERT_TRUE(data1.Empty());
        data1.AddPhoton(-0.3, TVector3(0.0, 0.0, -1.0), 1);
        ASSERT_TRUE(data1.Empty());

        PhotonCount data2 = CopyEmpty();
        ASSERT_TRUE(data2.Empty());
        data2.AddPhoton(0.2, TVector3(0.0, 1.0, 0.0), 1);
        ASSERT_TRUE(data2.Empty());

        PhotonCount data3 = CopyEmpty();
        ASSERT_TRUE(data3.Empty());
        PhotonCount::Iterator iter = data2.GetIterator();
        iter.Next();
        data3.AddNoise(1e6, iter);
        ASSERT_FALSE(data3.Empty());
    }

    /*
     * Make sure the time of a particular bin is calculated correctly.
     */
    TEST_F(DataStructuresTest, TestBinTime)
    {
        PhotonCount data = CopyEmpty();
        ASSERT_EQ(1, data.Bin(0.15));
        ASSERT_EQ(7, data.Bin(0.75));
    }

    /*
     * An out_of_range exception should be thrown if an invalid time is passed to the Bin() function.
     */
    TEST_F(DataStructuresTest, TestBinOutOfRange)
    {
        PhotonCount data = CopyEmpty();
        try
        {
            data.Bin(-0.1);
            FAIL() << "Exception not thrown";
        }
        catch(out_of_range& err)
        {
            ASSERT_EQ(string("Invalid time"), err.what());
        }
        try
        {
            data.Bin(0.99);
            FAIL() << "Exception not thrown";
        }
        catch(out_of_range& err)
        {
            ASSERT_EQ(string("Invalid time"), err.what());
        }
    }

    /*
     * Make sure the bin of a particular time is calculated correctly.
     */
    TEST_F(DataStructuresTest, TestTimeBin)
    {
        PhotonCount data = CopyEmpty();
        ASSERT_TRUE(Helper::ValuesEqual(0.15, data.Time(1), 1e-6));
        ASSERT_TRUE(Helper::ValuesEqual(0.75, data.Time(7), 1e-6));
    }

    /*
     * An out_of_range exception should be thrown if an invalid bin is passed to the Time() function.
     */
    TEST_F(DataStructuresTest, TestTimeOutOfRange)
    {
        PhotonCount data = CopyEmpty();
        try
        {
            data.Time(10);
            FAIL() << "Exception not thrown";
        }
        catch(exception& err)
        {
            ASSERT_EQ(string("Invalid bin"), err.what());
        }
        try
        {
            data.Time(-1);
            FAIL() << "Exception not thrown";
        }
        catch(exception& err)
        {
            ASSERT_EQ(string("Invalid bin"), err.what());
        }
    }

    /*
     * Check that the detector axis angle is calculated correctly.
     */
    TEST_F(DataStructuresTest, DetectorAxisAngle)
    {
        PhotonCount data = CopyEmpty();
        ASSERT_EQ(0.16, data.DetectorAxisAngle());
    }

    /*
     * Check that the pixel direction is calculated correctly using the iterator. The unit position of a pixel should be
     * its negative direction.
     */
    TEST_F(DataStructuresTest, PixelDirection)
    {
        PhotonCount data = CopyEmpty();
        PhotonCount::Iterator iter = data.GetIterator();
        iter.Next();
        iter.Next();
        iter.Next();
        iter.Next();

        TVector3 direction = TVector3(0, 0, 1);
        direction.RotateX(0.04);
        direction.RotateY(-0.04);
        ASSERT_TRUE(Helper::VectorsEqual(direction, data.Direction(iter), 1e-3));

        data.AddPhoton(0.45, -data.Direction(iter), 1);
        ASSERT_EQ(1, data.SumBins(iter));
    }

    /*
     * Check that the correct signal time series is returned from Signal()
     */
    TEST_F(DataStructuresTest, PixelSignal)
    {
        PhotonCount data = CopySample();
        PhotonCount::Iterator iter = data.GetIterator();

        iter.Next();
        iter.Next();
        iter.Next();
        ASSERT_EQ(Short1D(10, 0), data.Signal(iter));

        iter.Next();
        Short1D expected = Short1D(10, 0);
        expected[3] = 5;
        expected[4] = 1;
        expected[9] = 8;
        ASSERT_EQ(expected, data.Signal(iter));
    }

    /*
     * Check that bins are summed correctly in SumBins().
     */
    TEST_F(DataStructuresTest, SumBins)
    {
        PhotonCount data = CopySample();

        PhotonCount::Iterator iter = data.GetIterator();
        iter.Next();
        iter.Next();
        iter.Next();
        ASSERT_EQ(0, data.SumBins(iter));

        iter.Next();
        ASSERT_EQ(14, data.SumBins(iter));
    }

    /*
     * Check that the filtered bin sum is correct from SumBinsFiltered().
     */
    TEST_F(DataStructuresTest, SumBinsFiltered)
    {
        PhotonCount data = CopySample();
        Bool3D mask = data.GetFalseMatrix();
        mask[1][1][3] = true;
        mask[1][1][4] = true;
        PhotonCount::Iterator iter = data.GetIterator();

        iter.Next();
        iter.Next();
        ASSERT_EQ(0, data.SumBinsFiltered(iter, mask));

        iter.Next();
        iter.Next();
        ASSERT_EQ(6, data.SumBinsFiltered(iter, mask));
    }

    /*
     * Check the AverageTime function. The function should throw a invalid_argument error if the channel is empty
     */
    TEST_F(DataStructuresTest, AverageTime)
    {
        PhotonCount data = CopySample();
        PhotonCount::Iterator iter = data.GetIterator();
        iter.Next();
        iter.Next();
        ASSERT_TRUE(Helper::ValuesEqual(0.75, data.AverageTime(iter), 1e-6));

        iter.Next();
        try
        {
            data.AverageTime(iter);
            FAIL() << "Exception not thrown";
        }
        catch(invalid_argument& err)
        {
            ASSERT_EQ(string("Channel is empty, division by zero"), err.what());
        }

        iter.Next();
        ASSERT_TRUE(Helper::ValuesEqual(0.7, data.AverageTime(iter), 1e-6));
    }

    /*
     * Make sure the error in the mean of a pixel is calculated correctly.
     */
    TEST_F(DataStructuresTest, TimeError)
    {
        PhotonCount data = CopySample();
        PhotonCount::Iterator iter = data.GetIterator();
        iter.Next();
        iter.Next();
        ASSERT_EQ(0.1 / Sqrt(12.0), data.TimeError(iter));

        iter.Next();
        try
        {
            data.TimeError(iter);
            FAIL() << "Exception not thrown";
        }
        catch(invalid_argument& err)
        {
            ASSERT_EQ(string("Channel is empty, division by zero"), err.what());
        }

        iter.Next();
        double mean = data.AverageTime(iter);
        double vari = (5.0 * Sq(0.35 - mean) + Sq(0.45 - mean) + 8.0 * Sq(0.95 - mean)) / 14.0;
        vari += Sq(0.1) / 12.0;
        ASSERT_TRUE(Helper::ValuesEqual(Sqrt(vari / 14.0), data.TimeError(iter), 1e-6));
    }

    /*
     * Test the GetFalseMatrix() function.
     */
    TEST_F(DataStructuresTest, GetFalseMatrix)
    {
        PhotonCount data = CopyEmpty();
        Bool3D matrix = data.GetFalseMatrix();
        for (int i = 0; i < data.Size(); i++)
            for (int j = 0; j < data.Size(); j++)
                for (int k = 0; k < data.NBins(); k++)
                    ASSERT_FALSE(matrix[i][j][k]);
    }

    /*
     * Test the AddPhoton method when the photon is valid.
     */
    TEST_F(DataStructuresTest, AddValidPhoton)
    {
        PhotonCount data = CopyEmpty();
        TVector3 direction = TVector3(0, 0, 1);
        direction.RotateX(0.12);
        direction.RotateY(-0.04);
        data.AddPhoton(0.45, -direction, 3);
        ASSERT_FALSE(data.Empty());
        PhotonCount::Iterator iter = data.GetIterator();
        iter.Next();
        iter.Next();
        iter.Next();
        ASSERT_EQ(3, data.SumBins(iter));
    }

    /*
     * See what happens if a photon with an invalid position is added. This should leave the underlying container
     * unmodified.
     */
    TEST_F(DataStructuresTest, AddInvalidPosition)
    {
        PhotonCount data = CopyEmpty();
        data.AddPhoton(0.45, TVector3(0.0, 1.0, 0.0), 1);
        ASSERT_TRUE(data.Empty());

        PhotonCount::Iterator iter = data.GetIterator();
        while (iter.Next())
            ASSERT_EQ(0, data.SumBins(iter));
    }

    /*
     * See what happens if a photon with an undefined (0, 0, 0) position is added. This should throw a
     * invalid_argument. Also make sure we're able to add photons at angle pi/2.
     */
    TEST_F(DataStructuresTest, AddUndefinedPosition)
    {
        PhotonCount data = CopyEmpty();
        try
        {
            data.AddPhoton(0.35, TVector3(0, 0, 0), 1);
            FAIL() << "Exception not thrown";
        }
        catch (invalid_argument& err)
        {
            ASSERT_EQ(string("Direction cannot be a zero vector"), err.what());
        }
        data.AddPhoton(0.35, TVector3(1, 1, 0), 1);
    }

    /*
     * See what happens if a photon with an invalid time is added.
     */
    TEST_F(DataStructuresTest, AddInvalidTime)
    {
        PhotonCount data = CopyEmpty();
        data.AddPhoton(-0.1, TVector3(0, 0, 1), 1);
        ASSERT_TRUE(data.Empty());

        PhotonCount::Iterator iter = data.GetIterator();
        while (iter.Next())
            ASSERT_EQ(0, data.SumBins(iter));
    }

    /*
     * Make sure the noise subtraction works correctly.
     */
    TEST_F(DataStructuresTest, SubtractNoise)
    {
        PhotonCount data = CopySample();
        double universal_rate = 1e4;
        auto rate = (int) FriendRealNoiseRate(data, universal_rate);
        if (rate < 1) FAIL() << "Test is only trivially passing, increase the universal rate";

        PhotonCount::Iterator iter = data.GetIterator();
        iter.Next();
        iter.Next();
        iter.Next();
        iter.Next();
        int expected = data.SumBins(iter) - rate * 10;
        data.Subtract(universal_rate, iter);
        ASSERT_EQ(expected, data.SumBins(iter));
    }

    /*
     * Tests the AboveThreshold function.
     */
    TEST_F(DataStructuresTest, AboveThreshold)
    {
        PhotonCount data = CopySample();
        PhotonCount::Iterator iter = data.GetIterator();

        iter.Next();
        ASSERT_EQ(Bool1D(10, true), data.AboveThreshold(iter, -1));
        ASSERT_EQ(Bool1D(10, false), data.AboveThreshold(iter, 0));

        iter.Next();
        iter.Next();
        iter.Next();
        Bool1D expected = Bool1D(10, false);
        expected[9] = true;
        ASSERT_EQ(expected, data.AboveThreshold(iter, 6));
    }

    /*
     * Make sure FindThreshold correctly sets the threshold based on Gaussian probabilities. Note that, for the sample
     * data set, a universal noise rate of 1e4 implies a real noise rate of 6.4. The integral of a Gaussian above three
     * sigma is 0.001349.
     */
    TEST_F(DataStructuresTest, FindThreshold)
    {
        PhotonCount data = CopySample();
        ASSERT_EQ(15, data.FindThreshold(1e4, 3));
    }

    /*
     * Test the Subset function.
     */
    TEST_F(DataStructuresTest, Subset)
    {
        PhotonCount data = CopySample();
        Bool3D mat = data.GetFalseMatrix();
        mat[1][1][3] = true;
        data.Subset(mat);

        PhotonCount::Iterator iter = data.GetIterator();
        while (iter.Next())
        {
            Short1D signal = data.Signal(iter);
            if (iter.X() == 1 && iter.Y() == 1)
            {
                ASSERT_EQ(5, signal[3]);
                ASSERT_EQ(0, signal[4]);
                ASSERT_EQ(0, signal[9]);
            }
            else
            {
                ASSERT_EQ(Short1D(10, 0), signal);
            }
        }
    }

    /*
     * Checks that the Trim() function correctly adjusts the min/max times and reduces the size of all arrays.
     */
    TEST_F(DataStructuresTest, Trim)
    {
        PhotonCount data = CopySample();
        data.Trim();
        ASSERT_EQ(7, data.NBins());
        ASSERT_TRUE(Helper::ValuesEqual(0.35, data.Time(0), 1e-6));
        ASSERT_TRUE(Helper::ValuesEqual(0.95, data.Time(6), 1e-6));
    }
}

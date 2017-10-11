// DataStructuresTest.cpp
//
// Author: Matthew Dutson
//
// Tests of DataStructures.h

#include <gtest/gtest.h>
#include <TFile.h>

#include "DataStructures.h"
#include "Analysis.h"

namespace cherenkov_simulator
{
    /*
     * Note: this class will be able to access private member of the PhotonCount and Iterator classes
     */
    class DataStructuresTest : public ::testing::Test
    {
        PhotonCount::Params test_params;
        PhotonCount empty_data;
        PhotonCount sample_data;

        virtual void SetUp()
        {
            test_params = PhotonCount::Params();
            test_params.n_pixels = 4;
            test_params.max_bytes = 4000000;
            test_params.bin_size = 0.1;
            test_params.angular_size = 0.08;
            test_params.linear_size = 2.5;

            empty_data = PhotonCount(test_params, 0.0, 1.0);

            sample_data = PhotonCount(test_params, 0.0, 1.0);
            sample_data.IncrementCell(1, 0, 2, 7);
            sample_data.IncrementCell(2, 1, 1, 3);
            sample_data.IncrementCell(3, 1, 1, 3);
            sample_data.IncrementCell(1, 1, 1, 4);
            sample_data.IncrementCell(8, 1, 1, 9);
        }

        virtual void TearDown()
        {

        }

    // TODO: Shouldn't test fixtures have access to these methods without making them public?
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
        PhotonCount data = PhotonCount(CopyParams(), 0.0, 1.0);
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
            // TODO: Do we need to store this?
            PhotonCount data = PhotonCount(params, 0.0, 1.0);
            FAIL() << "Exception not thrown";
        }
        catch(std::exception& err)
        {
            ASSERT_EQ(std::string("Number of pixels must be even"), err.what());
        }
    }

    /*
     * Check what happens if we exceed the maximum number of bytes.
     */
    TEST_F(DataStructuresTest, OverMaxBytes)
    {
        PhotonCount::Params params = CopyParams();
        try
        {
            params.max_bytes = 10;
            PhotonCount data = PhotonCount(params, 0.0, 1.0);
            FAIL() << "Exception not thrown";
        }
        catch(std::exception& err)
        {
            ASSERT_EQ(std::string("Warning: too much memory requested due to shower direction"), err.what());
        }
        try
        {
            // TODO: May want to drop this to 340.
            params.max_bytes = 352;
            PhotonCount data = PhotonCount(params, 0.0, 1.0);
        }
        catch(std::exception& err)
        {
            FAIL() << "Exception thrown: " << err.what();
        }
    }

    /*
     * See what happens if the number of pixels passed to the constructor is zero. This shouldn't cause a division by
     * zero or weird behavior, so it should pass through.
     */
    TEST_F(DataStructuresTest, ZeroPixels)
    {
        PhotonCount::Params params = CopyParams();
        params.n_pixels = 0;
        try
        {
            PhotonCount data = PhotonCount(params, 0.0, 1.0);
        }
        catch(std::exception& err)
        {
            FAIL() << "Exception thrown: " << err.what();
        }
    }

    /*
     * See what happens if the bin size is zero. This will cause division by zero when determining the temporal bin, so
     * an exception should be thrown.
     */
    TEST_F(DataStructuresTest, ZeroBinSize)
    {
        PhotonCount::Params params = CopyParams();
        params.bin_size = 0;
        try
        {
            PhotonCount data = PhotonCount(params, 0.0, 1.0);
            FAIL() << "Exception not thrown";
        }
        catch(std::exception& err)
        {
            ASSERT_EQ(std::string("Bin size cannot be zero"), err.what());
        }
    }

    /*
     * See what happens if the angular size is zero (angular size is for a single pixel). This will cause division by
     * zero when determining the spatial bin, so an exception should be thrown.
     */
    TEST_F(DataStructuresTest, ZeroAngularSize)
    {
        PhotonCount::Params params = CopyParams();
        params.angular_size = 0;
        try
        {
            PhotonCount data = PhotonCount(params, 0.0, 1.0);
            FAIL() << "Exception not thrown";
        }
        catch(std::exception& err)
        {
            ASSERT_EQ(std::string("Angular size cannot be zero"), err.what());
        }
    }

    /*
     * See what happens if the linear size is zero (linear size is for a single pixel). Setting this to zero does not
     * necessarily cause division by zero, but does result in undefined behavior, so an exception should be thrown.
     */
    TEST_F(DataStructuresTest, ZeroLinearSize)
    {
        PhotonCount::Params params = CopyParams();
        params.linear_size = 0;
        try
        {
            // TODO: Do we need to actually declare a variable, or can we just call the constructor?
            PhotonCount data = PhotonCount(params, 0.0, 1.0);
            FAIL() << "Exception not thrown";
        }
        catch(std::exception& err)
        {
            ASSERT_EQ(std::string("Linear size cannot be zero"), err.what());
        }
    }

    /*
     * See what happens if we make the field of view very wide (close to 360 degrees).
     */
    TEST_F(DataStructuresTest, VeryWideField)
    {
        // TODO: I can see that this likely a problem, but not quite sure how to test it yet. Or maybe it's not a problem.
        FAIL() << "Not implemented";
    }

    /*
     * Checks that the correct pixels are valid.
     */
    TEST_F(DataStructuresTest, ValidPixels)
    {
        // TODO: May want to do a more complex example
        PhotonCount data = CopyEmpty();
        Bool2D valid = data.GetValid();
        for (int i = 0; i < valid.size(); i++)
        {
            for (int j = 0; j < valid.size(); j++)
            {
                if (TMath::Abs(i - 1.5) < 1 || TMath::Abs(j - 1.5) < 1)
                    ASSERT_TRUE(valid[i][j]);
                else
                    ASSERT_FALSE(valid[i][j]);
            }
        }
    }

    /*
     * Test the PhotonCount::Iterator. Steps through y and then steps through x. An exception should be thrown if the
     * iterator is at an invalid position and X() or Y() are called.
     */
    TEST_F(DataStructuresTest, TestIterator)
    {
        PhotonCount data = CopyEmpty();
        PhotonCount::Iterator iter = data.GetIterator();
        try
        {
            iter.X();
            FAIL() << "Exception not thrown";
        }
        catch(std::exception& err)
        {
            ASSERT_EQ(std::string("Call Next() before checking the iterator position"), err.what());
        }
        try
        {
            iter.Y();
            FAIL() << "Exception not thrown";
        }
        catch(std::exception& err)
        {
            ASSERT_EQ(std::string("Call Next() before checking the iterator position"), err.what());
        }

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
        PhotonCount data = CopyEmpty();
        PhotonCount::Iterator iter = data.GetIterator();

        iter.Next();
        iter.Next();
        iter.Reset();

        try
        {
            iter.X();
            FAIL() << "Exception not thrown";
        }
        catch(std::exception& err)
        {
            ASSERT_EQ(std::string("Call Next() before checking the iterator position"), err.what());
        }
        try
        {
            iter.Y();
            FAIL() << "Exception not thrown";
        }
        catch(std::exception& err)
        {
            ASSERT_EQ(std::string("Call Next() before checking the iterator position"), err.what());
        }

        ASSERT_TRUE(iter.Next());
        ASSERT_EQ(0, iter.X());
        ASSERT_EQ(1, iter.Y());

        ASSERT_TRUE(iter.Next());
        ASSERT_EQ(0, iter.X());
        ASSERT_EQ(2, iter.Y());
    }

    /*
     * Make sure the empty flag is updated correctly.
     */
    TEST_F(DataStructuresTest, TestEmpty)
    {
        PhotonCount data0 = CopyEmpty();
        ASSERT_TRUE(data0.Empty());
        data0.AddPhoton(0.2, TVector3(0.0, 0.0, 1.0), 1);
        ASSERT_FALSE(data0.Empty());

        // Emptiness shouldn't change if the photon is outside valid time bounds.
        PhotonCount data1 = CopyEmpty();
        ASSERT_TRUE(data1.Empty());
        data1.AddPhoton(-0.3, TVector3(0.0, 0.0, 1.0), 1);
        ASSERT_FALSE(data1.Empty());

        // Emptiness shouldn't change if the photon is outside valid spatial bounds.
        PhotonCount data2 = CopyEmpty();
        ASSERT_TRUE(data2.Empty());
        data2.AddPhoton(0.2, TVector3(0.0, 0.0, -1.0), 1);
        ASSERT_FALSE(data2.Empty());

        // Emptiness
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
     * Check what happens if the bin passed to Time(int) is out of range.
     */
    TEST_F(DataStructuresTest, TestBinOutOfRange)
    {
        PhotonCount data = CopyEmpty();
        try
        {
            data.Bin(-0.1);
            FAIL() << "Exception not thrown";
        }
        catch(std::exception& err)
        {
            ASSERT_EQ(std::string("Invalid time"), err.what());
        }
        try
        {
            data.Bin(1.1);
            FAIL() << "Exception not thrown";
        }
        catch(std::exception& err)
        {
            ASSERT_EQ(std::string("Invalid time"), err.what());
        }
    }

    /*
     * Make sure the bin of a particular time is calculated correctly.
     */
    TEST_F(DataStructuresTest, TestTimeBin)
    {
        PhotonCount data = CopyEmpty();
        ASSERT_EQ(0.1, data.Time(1));
        ASSERT_EQ(0.7, data.Time(7));
    }

    /*
     * Check what happens if the time passed to Bin(double) is out of range.
     */
    TEST_F(DataStructuresTest, TestTimeOutOfRange)
    {
        PhotonCount data = CopyEmpty();
        try
        {
            // TODO: Might want to change this to 9
            data.Time(10);
            FAIL() << "Exception not thrown";
        }
        catch(std::exception& err)
        {
            ASSERT_EQ(std::string("Invalid time"), err.what());
        }
        try
        {
            data.Time(-1);
            FAIL() << "Exception not thrown";
        }
        catch(std::exception& err)
        {
            ASSERT_EQ(std::string("Invalid time"), err.what());
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
     * Check that the pixel direction is calculated correctly using the iterator.
     */
    TEST_F(DataStructuresTest, PixelDirection)
    {
        FAIL() << "Not implemented";
    }

    /*
     * See what happens if an invalid iterator is passed to Direction. The Direction() method will check that they
     * both have the same size, but not that they both have the same validity mask, as this would be expensive.
     */
    TEST_F(DataStructuresTest, InvalidIterator)
    {
        PhotonCount::Params params = CopyParams();
        params.n_pixels = 2;
        PhotonCount data1 = PhotonCount(params, 0.0, 1.0);
        PhotonCount::Iterator iter1 = data1.GetIterator();
        iter1.Next();
        PhotonCount data2 = CopyEmpty();
        try
        {
            data2.Direction(iter1);
            FAIL() << "Exception not thrown";
        }
        catch(std::exception& err)
        {
            ASSERT_EQ(std::string("Iterator has invalid size"), err.what());
        }
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
        // TODO: May want to change this to 10.
        ASSERT_EQ(Short1D(11, 0), data.Signal(iter));

        iter.Next();
        // TODO: May want to change this to 10.
        Short1D expected = Short1D(11, 0);
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
        // TODO: May want to change this to 10
        mask[0][2] = Bool1D(11, false);
        mask[1][1][3] = true;
        mask[1][1][4] = true;
        mask[1][1][9] = false;
        PhotonCount::Iterator iter = data.GetIterator();
        iter.Next();
        iter.Next();
        ASSERT_EQ(0, data.SumBinsFiltered(iter, &mask));

        iter.Next();
        iter.Next();
        ASSERT_EQ(5, data.SumBinsFiltered(iter, &mask));
    }

    /*
     * Check the AverageTime function.
     */
    TEST_F(DataStructuresTest, AverageTime)
    {
        PhotonCount data = CopySample();
        PhotonCount::Iterator iter = data.GetIterator();
        iter.Next();
        iter.Next();
        ASSERT_EQ(0.7, data.AverageTime(iter));

        iter.Next();
        try
        {
            // TODO: Do we need to store this?
            double avg = data.AverageTime(iter);
            FAIL() << "Exception not thrown";
        }
        catch(std::exception& err)
        {
            ASSERT_EQ(std::string("Channel is empty, division by zero"), err.what());
        }

        iter.Next();
        ASSERT_EQ(0.65, data.AverageTime(iter));
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
        ASSERT_EQ(0.1 / TMath::Sqrt(12.0), data.TimeError(iter));

        iter.Next();
        try
        {
            // TODO: Do we need to store this?
            double var = data.TimeError(iter);
            FAIL() << "Exception not thrown";
        }
        catch(std::exception& err)
        {
            ASSERT_EQ(std::string("Channel is empty, division by zero"), err.what());
        }

        // TODO: Is our data centrally binned or left binned? This will likely impact our estimate of errors.
        iter.Next();
        double expected = TMath::Sqrt(117.5 / 14.0 + TMath::Sq(0.1) / 12.0);
        ASSERT_EQ(expected, data.TimeError(iter));
    }

    /*
     * Test the GetFalseMatrix() function.
     */
    TEST_F(DataStructuresTest, GetFalseMatrix)
    {
        PhotonCount data = CopyEmpty();
        Bool3D matrix = data.GetFalseMatrix();
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                // TODO: May need to change this to 10
                for (int k = 0; k < 11; k++)
                {
                    ASSERT_FALSE(matrix[i][j][k]);
                }
            }
        }
    }

    /*
     * Test the AddPhoton method when the photon is valid.
     */
    TEST_F(DataStructuresTest, AddValidPhoton)
    {
        FAIL() << "Not implemented";
    }

    /*
     * See what happens if a photon with an invalid position is added. This should leave the underlying container
     * unmodified.
     */
    TEST_F(DataStructuresTest, AddInvalidPosition)
    {
        PhotonCount data = CopyEmpty();
        data.AddPhoton(0.35, TVector3(0, 0, -1), 1);
        ASSERT_TRUE(data.Empty());

        PhotonCount::Iterator iter = data.GetIterator();
        while (iter.Next())
            ASSERT_EQ(0, data.SumBins(iter));
    }

    /*
     * See what happens if a photon with an undefined (0, 0, 0) position is added. This should throw a
     * std::invalid_argument. Also make sure we're able to add photons at angle pi/2.
     */
    TEST_F(DataStructuresTest, AddUndefinedPosition)
    {
        PhotonCount data = CopyEmpty();
        try
        {
            data.AddPhoton(0.35, TVector3(0, 0, 0), 1);
            FAIL() << "Exception not thrown";
        }
        catch (std::exception& err)
        {
            ASSERT_EQ(std::string("Direction cannot be a zero vector"), err.what());
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
     * Make sure the noise addition does basically what we would expect (this is a probabilistic test).
     */
    TEST_F(DataStructuresTest, AddNoise)
    {
        FAIL() << "Not implemented";
    }

    /*
     * Make sure the noise subtraction works correctly.
     */
    TEST_F(DataStructuresTest, SubtractNoise)
    {
        FAIL() << "Not implemented";
    }

    /*
     * Tests the AboveThreshold function.
     */
    TEST_F(DataStructuresTest, AboveThreshold)
    {
        FAIL() << "Not implemented";
    }

    /*
     * Make sure FindThreshold correctly sets the threshold based on Gaussian probabilities.
     */
    TEST_F(DataStructuresTest, FindThreshold)
    {
        FAIL() << "Not implemented";
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
                // TODO: May want to change this to 10
                ASSERT_EQ(Short1D(11, 0), signal);
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
        ASSERT_EQ(0.3, data.Time(0));
        ASSERT_EQ(0.9, data.Time(6));
    }
}

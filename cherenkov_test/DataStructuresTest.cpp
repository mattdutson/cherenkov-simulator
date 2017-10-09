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
        PhotonCount test_data;

        virtual void SetUp()
        {
            test_params = PhotonCount::Params();
            test_params.n_pixels = 4;
            test_params.max_bytes = 4000000;
            test_params.bin_size = 0.1;
            test_params.angular_size = 0.08;
            test_params.linear_size = 2.5;

            test_data = PhotonCount(test_params, 0.0, 1.0);
        }

        virtual void TearDown()
        {

        }

    // TODO: Shouldn't test fixtures have access to these methods without making them public?
    public:

        PhotonCount CopyData()
        {
            return test_data;
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
        PhotonCount data = CopyData();
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
        PhotonCount data = CopyData();
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
        PhotonCount data = CopyData();
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
        PhotonCount data0 = CopyData();
        ASSERT_TRUE(data0.Empty());
        data0.AddPhoton(0.2, TVector3(0.0, 0.0, 1.0), 1);
        ASSERT_FALSE(data0.Empty());

        // Emptiness shouldn't change if the photon is outside valid time bounds.
        PhotonCount data1 = CopyData();
        ASSERT_TRUE(data1.Empty());
        data1.AddPhoton(-0.3, TVector3(0.0, 0.0, 1.0), 1);
        ASSERT_FALSE(data1.Empty());

        // Emptiness shouldn't change if the photon is outside valid spatial bounds.
        PhotonCount data2 = CopyData();
        ASSERT_TRUE(data2.Empty());
        data2.AddPhoton(0.2, TVector3(0.0, 0.0, -1.0), 1);
        ASSERT_FALSE(data2.Empty());

        // Emptiness
        PhotonCount data3 = CopyData();
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
        PhotonCount data = CopyData();
        ASSERT_EQ(1, data.Bin(0.15));
        ASSERT_EQ(7, data.Bin(0.75));
    }

    /*
     * Check what happens if the bin passed to Time(int) is out of range.
     */
    TEST_F(DataStructuresTest, TestBinOutOfRange)
    {
        PhotonCount data = CopyData();
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
        PhotonCount data = CopyData();
        ASSERT_EQ(0.1, data.Time(1));
        ASSERT_EQ(0.7, data.Time(7));
    }

    /*
     * Check what happens if the time passed to Bin(double) is out of range.
     */
    TEST_F(DataStructuresTest, TestTimeOutOfRange)
    {
        PhotonCount data = CopyData();
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
        PhotonCount data = CopyData();
        ASSERT_EQ(0.16, data.DetectorAxisAngle());
    }

    /*
     * Check that the pixel direction is calculated correctly using the iterator.
     */
    TEST_F(DataStructuresTest, PixelDirection)
    {

    }

    /*
     * See what happens if an invalid iterator is passed to Direction.
     */
    TEST_F(DataStructuresTest, InvalidIterator)
    {

    }

    /*
     * Check that the correct signal time series is returned from Signal()
     */
    TEST_F(DataStructuresTest, PixelSignal)
    {

    }

    /*
     * Check that bins are summed correctly in SumBins().
     */
    TEST_F(DataStructuresTest, SumBins)
    {

    }

    /*
     * Check that the filtered bin sum is correct from SumBinsFiltered().
     */
    TEST_F(DataStructuresTest, SumBinsFiltered)
    {

    }

    /*
     * Check the AverageTime function.
     */
    TEST_F(DataStructuresTest, AverageTime)
    {

    }

    /*
     * Make sure the error in the mean of a pixel is calculated correctly.
     */
    TEST_F(DataStructuresTest, TimeError)
    {

    }

    /*
     * Test the GetFalseMatrix() function.
     */
    TEST_F(DataStructuresTest, GetFalseMatrix)
    {

    }

    /*
     * Test the AddPhoton method when the photon is valid.
     */
    TEST_F(DataStructuresTest, AddValidPhoton)
    {

    }

    /*
     * See what happens if a photon with an invalid time is added.
     */
    TEST_F(DataStructuresTest, AddInvalidPosition)
    {

    }

    /*
     * See what happens if a photon with an invalid position is added.
     */
    TEST_F(DataStructuresTest, AddInvalidTime)
    {

    }

    /*
     * Make sure the noise addition does basically what we would expect (this is a probabilistic test).
     */
    TEST_F(DataStructuresTest, AddNoise)
    {

    }

    /*
     * Make sure the noise subtraction works correctly.
     */
    TEST_F(DataStructuresTest, SubtractNoise)
    {

    }

    /*
     * Tests the AboeThreshold function.
     */
    TEST_F(DataStructuresTest, AboveThreshold)
    {

    }

    /*
     * Make sure FindThreshold correctly sets the threshold based on Gaussian probabilities.
     */
    TEST_F(DataStructuresTest, FindThreshold)
    {

    }

    /*
     * Test the Subset function.
     */
    TEST_F(DataStructuresTest, Subset)
    {

    }

    /*
     * Checks that the Trim() function correctly adjusts the min/max times and reduces the size of all arrays.
     */
    TEST_F(DataStructuresTest, Trim)
    {

    }
}

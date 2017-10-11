// GeometricTest.cpp
//
// Author: Matthew Dutson
//
// Tests of Geometric.h

#include <gtest/gtest.h>
#include <TMath.h>

#include "Geometric.h"
#include "Helper.h"

using namespace TMath;

namespace cherenkov_simulator
{

/*
 * Note: this class will be able to access private member of the Plane, Shower, and Ray classes.
 */
class GeometricTest : public ::testing::Test
{
private:
    
    Ray test_ray_1;
    Ray test_ray_2;

    virtual void SetUp()
    {
        test_ray_1 = Ray(TVector3(7, 2, 3), TVector3(8, 8, 8), 1.7);
        test_ray_2 = Ray(TVector3(0, 0, 0), TVector3(-1, 2, -1), -0.8);
    }

    virtual void TearDown()
    {
    }
    
public:
    
    Ray CopyRay1()
    {
        return test_ray_1;
    }

    Ray CopyRay2()
    {
        return test_ray_2;
    }
};
    /*
     * Test the default constructor for a plane object.
     */
    TEST_F(GeometricTest, DefaultPlane)
    {
        Plane plane = Plane();
        ASSERT_EQ(plane.Normal(), TVector3(0, 0, 1));
        ASSERT_EQ(0, plane.Coefficient());
    }

    /*
     * Check that the normal vector is unit.
     */
    TEST_F(GeometricTest, PlaneUnitNormal)
    {
        Plane plane = Plane(TVector3(0, 0, 12), TVector3(0, 0, 0));
        ASSERT_EQ(TVector3(0, 0, 1), plane.Normal());
    }

    /*
     * Test that the correct value of the coefficient is returned.
     */
    TEST_F(GeometricTest, PlaneCoefficient)
    {
        Plane plane = Plane(TVector3(0, 0, 1), TVector3(45, -12, 33));
        ASSERT_EQ(33, plane.Coefficient());
    }

    /*
     * Set a non-default normal vector through the constructor.
     */
    TEST_F(GeometricTest, SetPlaneNormal)
    {
        Plane plane = Plane(TVector3(1, 1, 1), TVector3(45, -12, 33));
        ASSERT_EQ(plane.Normal(), TVector3(1 / Sqrt(3), 1 / Sqrt(3), 1 / Sqrt(3)));
        ASSERT_EQ(66, plane.Coefficient());
    }

    /*
     * Ensure that an exception is thrown if an invalid direction vector is passed to the constructor.
     */
    TEST_F(GeometricTest, BadNormalDirection)
    {
        try
        {
            Ray ray = Ray(TVector3(10, 120, -11), TVector3(0, 0, 0), 0.04);
            FAIL() << "Exception not thrown";
        }
        catch(std::exception& err)
        {
            ASSERT_EQ(std::string("Plane normal vector must be nonzero"), err.what());
        }
    }

    /*
     * Check that the InFrontOf function operates correctly. Should return false if the plane is fixed at the origin.
     */
    TEST_F(GeometricTest, InFrontOf)
    {
        Plane front = Plane(TVector3(1, 0, 0), TVector3(10, 0, 0));
        Plane origin = Plane(TVector3(1, 0, 0), TVector3(0, 0, 0));
        Plane back = Plane(TVector3(1, 0, 0), TVector3(-10, 0, 0));
        ASSERT_TRUE(front.InFrontOf(TVector3(1, 0, 0)));
        ASSERT_FALSE(origin.InFrontOf(TVector3(1, 0, 0)));
        ASSERT_FALSE(back.InFrontOf(TVector3(1, 0, 0)));
    }

    /*
     * Tests the default constructor for the Ray class.
     */
    TEST_F(GeometricTest, DefaultRay)
    {
        Ray ray = Ray();
        ASSERT_EQ(TVector3(0, 0, 1), ray.Direction());
        ASSERT_EQ(TVector3(0, 0, 0), ray.Position());
        ASSERT_EQ(0, ray.Time());
    }

    /*
     * Tests the non-default constructor for the Ray class.
     */
    TEST_F(GeometricTest, UserRay)
    {
        Ray ray = CopyRay1();
        ASSERT_EQ(TVector3(7, 2, 3), ray.Position());
        ASSERT_EQ(TVector3(1, 1, 1).Unit(), ray.Direction());
        ASSERT_EQ(1.7, ray.Time());
    }

    /*
     * Tests the behavior when an invalid direction (0, 0, 0) is passed to the Ray constructor.
     */
    TEST_F(GeometricTest, InvalidDirection)
    {
        try
        {
            // TODO: Do we need to set this?
            Ray ray = Ray(TVector3(0, 0, 0), TVector3(0, 0, 0), 0);
            FAIL() << "Exception not thrown";
        }
        catch(std::exception& err)
        {
            ASSERT_EQ(std::string("Ray direction must be nonzero"), err.what());
        }
    }

    /*
     * Tests that the Ray's position is returned correctly and updated when a propagate function is called.
     */
    TEST_F(GeometricTest, TestPosition)
    {
        Ray ray = CopyRay1();
        ASSERT_EQ(TVector3(7, 2, 3), ray.Position());

        Plane plane = Plane(TVector3(1, 0, 0), TVector3(10, 0, 0));
        ray.PropagateToPlane(plane);
        ASSERT_EQ(TVector3(10, 5, 6), ray.Position());

        ray.PropagateToPoint(TVector3(12, 7, 3));
        ASSERT_EQ(TVector3(12, 7, 3), ray.Position());
    }

    /*
     * Checks that the velocity is given the correct magnitude (speed of light), and that it is updated on a direction
     * change.
     */
    TEST_F(GeometricTest, TestVelocity)
    {
        Ray ray = CopyRay1();
        ASSERT_EQ(TVector3(1, 1, 1).Unit() * Utility::c_cent, ray.Velocity());

        ray.Reflect(TVector3(-1, 0, 0));
        ASSERT_EQ(TVector3(-1, 1, 1).Unit() * Utility::c_cent, ray.Velocity());

        ray.SetDirection(TVector3(-1, 2, 0));
        ASSERT_EQ(TVector3(-1, 2, 0).Unit() * Utility::c_cent, ray.Velocity());
    }

    /*
     * Checks that the direction has unit magnitude and is updated on a direction change.
     */
    TEST_F(GeometricTest, TestDirection)
    {
        Ray ray = CopyRay1();
        ASSERT_EQ(TVector3(1, 1, 1).Unit(), ray.Direction());

        ray.Reflect(TVector3(-1, 0, 0));
        ASSERT_EQ(TVector3(-1, 1, 1).Unit(), ray.Direction());

        ray.SetDirection(TVector3(-1, 2, 0));
        ASSERT_EQ(TVector3(-1, 2, 0).Unit(), ray.Direction());
    }

    /*
     * Checks that an exception is thrown if the SetDirection function is passed a zero vector.
     */
    TEST_F(GeometricTest, SetInvalidDirection)
    {
        Ray ray = CopyRay1();
        try
        {
            ray.SetDirection(TVector3(0, 0, 0));
            FAIL() << "Exception not thrown";
        }
        catch(std::exception& err)
        {
            ASSERT_EQ(std::string("Direction vector must be nonzero"), err.what());
        }
    }

    /*
     * Checks that the time is returned correctly, and is updated when the Ray propagates.
     */
    TEST_F(GeometricTest, Time)
    {
        Ray ray1 = CopyRay1();
        ASSERT_EQ(1.7, ray1.Time());

        Plane plane = Plane(TVector3(1, 0, 0), TVector3(10, 0, 0));
        ray1.PropagateToPlane(plane);
        ASSERT_EQ(1.7 + TVector3(3, 3, 3).Mag() / Utility::c_cent, ray1.Time());

        Ray ray2 = CopyRay1();
        ray2.PropagateToPoint(TVector3(-1, 7, 8));
        ASSERT_EQ(1.7 + TVector3(-8, 5, 5).Mag() / Utility::c_cent, ray2.Time());
    }

    /*
     * Tests the PropagateToPoint method when the point is along the direction of the Ray.
     */
    TEST_F(GeometricTest, PropagateToPoint)
    {
        Ray ray = CopyRay1();
        TVector3 dir_init = ray.Direction();
        ray.PropagateToPoint(TVector3(20, 15, 16));
        ASSERT_EQ(TVector3(20, 15, 16), ray.Position());
        ASSERT_EQ(dir_init, ray.Direction());
        ASSERT_EQ(1.7 + TVector3(13, 13, 13).Mag() / Utility::c_cent, ray.Time());
    }

    /*
     * Tests the PropagateToPoint method when the point is not along the direction of the Ray (the direction should be
     * changed and it should then be propagated).
     */
    TEST_F(GeometricTest, PropagateToPointChange)
    {
        Ray ray = CopyRay1();
        ray.PropagateToPoint(TVector3(-8, 97, 4));
        ASSERT_EQ(TVector3(-8, 97, 4), ray.Position());
        ASSERT_EQ(TVector3(-15, 55, 1).Unit(), ray.Direction());
        ASSERT_EQ(1.7 + TVector3(-15, 55, 1).Mag() / Utility::c_cent, ray.Time());
    }

    /*
     * Tests the PropagateToPlane method when the plane is in front of or behind the Ray.
     */
    TEST_F(GeometricTest, PropagateToPlane)
    {
        Ray ray = CopyRay2();
        
    }

    /*
     * Tests the PropagateToPlane method when the plane is parallel to the Ray.
     */
    TEST_F(GeometricTest, PropagateToPlaneParallel)
    {

    }

    /*
     * Tests the PlaneImpact method when the plane is in front of or behind the Ray.
     */
    TEST_F(GeometricTest, PlaneImpact)
    {

    }

    /*
     * Tests the PlaneImpact method when the plane is parallel to the Ray (the Ray's current position should be
     * returned).
     */
    TEST_F(GeometricTest, PlaneImpactParallel)
    {

    }

    /*
     * Test the TimeToPlane method when the plane is in front of or behind the Ray.
     */
    TEST_F(GeometricTest, TimeToPlane)
    {

    }

    /*
     * Tests the TimeToPlane method when the plane is parallel to the Ray (infinity should be returned).
     */
    TEST_F(GeometricTest, TimeToPlaneParallel)
    {

    }

    /*
     * Tests the Reflect method. The sign of the normal vector SHOULD NOT matter.
     */
    TEST_F(GeometricTest, Reflect)
    {

    }

    /*
     * Tests the Reflect method when a zero vector is passed as the normal. In this case, std::invalid_argument should
     * be thrown.
     */
    TEST_F(GeometricTest, ReflectZeroNormal)
    {

    }

    /*
     * Tests a normal application of the Refract method, both from less dense to more dense and vice versa.
     */
    TEST_F(GeometricTest, Refract)
    {

    }

    /*
     * Tests the Refract method when a zero vector is passed as the normal. In this case, std::invalid_argument should
     * be thrown.
     */
    TEST_F(GeometricTest, RefractZeroNormal)
    {

    }

    /*
     * Tests the Refract method when total internal reflection occurs (will only be the case if n_out < n_in).
     */
    TEST_F(GeometricTest, RefractCritAngle)
    {

    }

    /*
     * Test the Refract method when an n < 1 is passed. In this case, std::invalid_argument should be thrown.
     */
    TEST_F(GeometricTest, RefractInvalidN)
    {

    }

    /*
     * Test the Transform method.
     */
    TEST_F(GeometricTest, Transform)
    {

    }

    /*
     * Test the default constructor for Shower.
     */
    TEST_F(GeometricTest, DefaultShower)
    {

    }

    /*
     * Test the non-default constructor for Shower.
     */
    TEST_F(GeometricTest, UserShower)
    {

    }

    /*
     * Test the Shower constructor when a non-positive energy is passed. In this case, std::invalid_argument should be
     * thrown.
     */
    TEST_F(GeometricTest, NonPositiveEnergy)
    {

    }

    /*
     * Test the Shower constructor when a non-positive XMax is passed. In this case, std::invalid_argument should be
     * thrown.
     */
    TEST_F(GeometricTest, NonPositiveXMax)
    {

    }

    /*
     * Test the Shower constructor when a non-positive NMax is passed. In this case, std::invalid_argument should be
     * thrown.
     */
    TEST_F(GeometricTest, NonPositiveNMax)
    {

    }

    /*
     * Test the Shower constructor when a non-positive Rho0 is passed. In this case, std::invalid_argument should be
     * thrown.
     */
    TEST_F(GeometricTest, NonPositiveRho0)
    {

    }

    /*
     * Test the Shower constructor when a non-positive scale height is passed. In this case, std::invalid_argument
     * should be thrown.
     */
    TEST_F(GeometricTest, NonPositiveScaleH)
    {

    }

    /*
     * Test the Shower constructor when a non-positive Delta0 is passed. In this case, std::invalid_argument should be
     * thrown.
     */
    TEST_F(GeometricTest, NonPositiveDelta0)
    {

    }

    /*
     * Test the Age function.
     */
    TEST_F(GeometricTest, ShowerAge)
    {

    }

    /*
     * Test the EnergyMeV function.
     */
    TEST_F(GeometricTest, EnergyMeV)
    {

    }

    /*
     * Test the EnergyeV function.
     */
    TEST_F(GeometricTest, EnergyeV)
    {

    }

    /*
     * Test the ImpactParam function.
     */
    TEST_F(GeometricTest, ImpactParam)
    {

    }

    /*
     * Test the ImpactAngle function.
     */
    TEST_F(GeometricTest, ImpactAngle)
    {

    }

    /*
     * Test the LocalRho function.
     */
    TEST_F(GeometricTest, LocalRho)
    {

    }

    /*
     * Test the LocalDelta function.
     */
    TEST_F(GeometricTest, LocalDelta)
    {

    }

    /*
     * Test the GaisserHillas function.
     */
    TEST_F(GeometricTest, GaisserHillas)
    {

    }

    /*
     * Test the EThresh function.
     */
    TEST_F(GeometricTest, EThresh)
    {

    }

    /*
     * Test the IncrementDepth function. Note that this is defined for zero and negative depths.
     */
    TEST_F(GeometricTest, IncrementDepth)
    {

    }
}
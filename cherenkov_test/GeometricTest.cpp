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
     * Test the default constructor for a plane object.
     */
    TEST(GeometricTest, DefaultPlane)
    {
        Plane plane = Plane();
        ASSERT_EQ(plane.Normal(), TVector3(0, 0, 1));
        ASSERT_EQ(0, plane.Coefficient());
    }

    /*
     * Check that the normal vector is unit.
     */
    TEST(GeometricTest, PlaneUnitNormal)
    {
        Plane plane = Plane(TVector3(0, 0, 12), TVector3(0, 0, 0));
        ASSERT_EQ(TVector3(0, 0, 1), plane.Normal());
    }

    /*
     * Test that the correct value of the coefficient is returned.
     */
    TEST(GeometricTest, PlaneCoefficient)
    {
        Plane plane = Plane(TVector3(0, 0, 1), TVector3(45, -12, 33));
        ASSERT_EQ(33, plane.Coefficient());
    }

    /*
     * Set a non-default normal vector through the constructor.
     */
    TEST(GeometricTest, SetPlaneNormal)
    {
        Plane plane = Plane(TVector3(1, 1, 1), TVector3(45, -12, 33));
        ASSERT_EQ(plane.Normal(), TVector3(1 / Sqrt(3), 1 / Sqrt(3), 1 / Sqrt(3)));
        ASSERT_EQ(66, plane.Coefficient());
    }

    /*
     * Ensure that an exception is thrown if an invalid direction vector is passed to the constructor.
     */
    TEST(GeometricTest, BadNormalDirection)
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
    TEST(GeometricTest, InFrontOf)
    {

    }

    /*
     * Tests the default constructor for the Ray class.
     */
    TEST(GeometricTest, DefaultRay)
    {

    }

    /*
     * Tests the non-default constructor for the Ray class.
     */
    TEST(GeometricTest, UserRay)
    {

    }

    /*
     * Tests the behavior when an invalid direction (0, 0, 0) is passed to the Ray constructor.
     */
    TEST(GeometricTest, InvalidDirection)
    {

    }

    /*
     * Tests that the Ray's position is returned correctly and updated when a propagate function is called.
     */
    TEST(GeometricTest, TestPosition)
    {

    }

    /*
     * Checks that the velocity is given the correct magnitude (speed of light), and that it is updated on a direction
     * change.
     */
    TEST(GeometricTest, TestVelocity)
    {

    }

    /*
     * Checks that the direction has unit magnitude and is updated on a direction change.
     */
    TEST(GeometricTest, TestDirection)
    {

    }

    /*
     * Checks that an exception is thrown if the SetDirection function is passed a zero vector.
     */
    TEST(GeometricTest, SetInvalidDirection)
    {

    }

    /*
     * Checks that the time is returned correctly, and is updated when the Ray propagates.
     */
    TEST(GeometricTest, Time)
    {

    }

    /*
     * Tests the PropagateToPoint method when the point is along the direction of the Ray.
     */
    TEST(GeometricTest, PropagateToPoint)
    {

    }

    /*
     * Tests the PropagateToPoint method when the point is not along the direction of the Ray (the direction should be
     * changed and it should then be propagated).
     */
    TEST(GeometricTest, PropagateToPointChange)
    {

    }

    /*
     * Tests the PropagateToPlane method when the plane is in front of or behind the Ray.
     */
    TEST(GeometricTest, PropagateToPlane)
    {

    }

    /*
     * Tests the PropagateToPlane method when the plane is parallel to the Ray.
     */
    TEST(GeometricTest, PropagateToPlaneParallel)
    {

    }

    /*
     * Tests the PlaneImpact method when the plane is in front of or behind the Ray.
     */
    TEST(GeometricTest, PlaneImpact)
    {

    }

    /*
     * Tests the PlaneImpact method when the plane is parallel to the Ray (the Ray's current position should be
     * returned).
     */
    TEST(GeometricTest, PlaneImpactParallel)
    {

    }

    /*
     * Test the TimeToPlane method when the plane is in front of or behind the Ray.
     */
    TEST(GeometricTest, TimeToPlane)
    {

    }

    /*
     * Tests the TimeToPlane method when the plane is parallel to the Ray (infinity should be returned).
     */
    TEST(GeometricTest, TimeToPlaneParallel)
    {

    }

    /*
     * Tests the Reflect method. The sign of the normal vector SHOULD NOT matter.
     */
    TEST(GeometricTest, Reflect)
    {

    }

    /*
     * Tests the Reflect method when a zero vector is passed as the normal. In this case, std::invalid_argument should
     * be thrown.
     */
    TEST(GeometricTest, ReflectZeroNormal)
    {

    }

    /*
     * Tests a normal application of the Refract method, both from less dense to more dense and vice versa.
     */
    TEST(GeometricTest, Refract)
    {

    }

    /*
     * Tests the Refract method when a zero vector is passed as the normal. In this case, std::invalid_argument should
     * be thrown.
     */
    TEST(GeometricTest, RefractZeroNormal)
    {

    }

    /*
     * Tests the Refract method when total internal reflection occurs (will only be the case if n_out < n_in).
     */
    TEST(GeometricTest, RefractCritAngle)
    {

    }

    /*
     * Test the Refract method when an n < 1 is passed. In this case, std::invalid_argument should be thrown.
     */
    TEST(GeometricTest, RefractInvalidN)
    {

    }

    /*
     * Test the Transform method.
     */
    TEST(GeometricTest, Transform)
    {

    }

    /*
     * Test the default constructor for Shower.
     */
    TEST(GeometricTest, DefaultShower)
    {

    }

    /*
     * Test the non-default constructor for Shower.
     */
    TEST(GeometricTest, UserShower)
    {

    }

    /*
     * Test the Shower constructor when a non-positive energy is passed. In this case, std::invalid_argument should be
     * thrown.
     */
    TEST(GeometricTest, NonPositiveEnergy)
    {

    }

    /*
     * Test the Shower constructor when a non-positive XMax is passed. In this case, std::invalid_argument should be
     * thrown.
     */
    TEST(GeometricTest, NonPositiveXMax)
    {

    }

    /*
     * Test the Shower constructor when a non-positive NMax is passed. In this case, std::invalid_argument should be
     * thrown.
     */
    TEST(GeometricTest, NonPositiveNMax)
    {

    }

    /*
     * Test the Shower constructor when a non-positive Rho0 is passed. In this case, std::invalid_argument should be
     * thrown.
     */
    TEST(GeometricTest, NonPositiveRho0)
    {

    }

    /*
     * Test the Shower constructor when a non-positive scale height is passed. In this case, std::invalid_argument
     * should be thrown.
     */
    TEST(GeometricTest, NonPositiveScaleH)
    {

    }

    /*
     * Test the Shower constructor when a non-positive Delta0 is passed. In this case, std::invalid_argument should be
     * thrown.
     */
    TEST(GeometricTest, NonPositiveDelta0)
    {

    }

    /*
     * Test the Age function.
     */
    TEST(GeometricTest, ShowerAge)
    {

    }

    /*
     * Test the EnergyMeV function.
     */
    TEST(GeometricTest, EnergyMeV)
    {

    }

    /*
     * Test the EnergyeV function.
     */
    TEST(GeometricTest, EnergyeV)
    {

    }

    /*
     * Test the ImpactParam function.
     */
    TEST(GeometricTest, ImpactParam)
    {

    }

    /*
     * Test the ImpactAngle function.
     */
    TEST(GeometricTest, ImpactAngle)
    {

    }

    /*
     * Test the LocalRho function.
     */
    TEST(GeometricTest, LocalRho)
    {

    }

    /*
     * Test the LocalDelta function.
     */
    TEST(GeometricTest, LocalDelta)
    {

    }

    /*
     * Test the GaisserHillas function.
     */
    TEST(GeometricTest, GaisserHillas)
    {

    }

    /*
     * Test the EThresh function.
     */
    TEST(GeometricTest, EThresh)
    {

    }

    /*
     * Test the IncrementDepth function. Note that this is defined for zero and negative depths.
     */
    TEST(GeometricTest, IncrementDepth)
    {

    }
}
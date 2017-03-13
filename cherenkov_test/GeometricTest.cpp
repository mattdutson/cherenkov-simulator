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
    TEST(GeometricTest, ConstructPlane)
    {
        Plane plane = Plane();
        ASSERT_EQ(plane.Normal(), TVector3(0, 0, 1));
        ASSERT_EQ(plane.Coefficient(), 0);
    }

    /*
     * Check that the normal vector is unit.
     */
    TEST(GeometricTest, PlaneUnitNormal)
    {
        Plane plane = Plane(TVector3(0, 0, 12), TVector3(0, 0, 0));
        ASSERT_EQ(plane.Normal(), TVector3(0, 0, 1));
    }

    /*
     * Test that the correct value of the coefficient is returned.
     */
    TEST(GeometricTest, PlaneCoefficient)
    {
        Plane plane = Plane(TVector3(0, 0, 1), TVector3(45, -12, 33));
        ASSERT_EQ(plane.Coefficient(), 33);
    }

    /*
     * Set a non-default normal vector through the constructor.
     */
    TEST(GeometricTest, SetPlaneNormal)
    {
        Plane plane = Plane(TVector3(1, 1, 1), TVector3(45, -12, 33));
        ASSERT_EQ(plane.Normal(), TVector3(1 / Sqrt(3), 1 / Sqrt(3), 1 / Sqrt(3)));
    }

    /*
     * Make sure the constructor corrects an invalid direction vector.
     */
    TEST(GeometricTest, BadNormalDirection)
    {
        Ray ray = Ray(TVector3(10, 120, -11), TVector3(0, 0, 0), 0.04);
        ASSERT_TRUE(Helper::VectorsEqual(ray.Direction(), TVector3(0, 0, 1), 0.001));
    }
}
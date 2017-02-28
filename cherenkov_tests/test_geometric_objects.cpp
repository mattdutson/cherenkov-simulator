//
// Created by Matthew Dutson on 1/17/17.
//

#include <gtest/gtest.h>
#include <TMath.h>

#include "geometric_objects.h"
#include "helper.h"

using namespace cherenkov_library;
using namespace TMath;

namespace cherenkov_tests
{
    /*
     * Test the default constructor for a plane object.
     */
    TEST(GeometricObjectsTest, ConstructPlane)
    {
        Plane plane = Plane();
        ASSERT_EQ(plane.Normal(), TVector3(0, 0, 1));
        ASSERT_EQ(plane.Coefficient(), 0);
    }

    /*
     * Check that the normal vector is unit.
     */
    TEST(GeometricObjectsTest, PlaneUnitNormal)
    {
        Plane plane = Plane(TVector3(0, 0, 12), TVector3(0, 0, 0));
        ASSERT_EQ(plane.Normal(), TVector3(0, 0, 1));
    }

    /*
     * Test that the correct value of the coefficient is returned.
     */
    TEST(GeometricObjectsTest, PlaneCoefficient)
    {
        Plane plane = Plane(TVector3(0, 0, 1), TVector3(45, -12, 33));
        ASSERT_EQ(plane.Coefficient(), 33);
    }

    /*
     * Set a non-default normal vector through the constructor.
     */
    TEST(GeometricObjectsTest, SetPlaneNormal)
    {
        Plane plane = Plane(TVector3(1, 1, 1), TVector3(45, -12, 33));
        ASSERT_EQ(plane.Normal(), TVector3(1 / Sqrt(3), 1 / Sqrt(3), 1 / Sqrt(3)));
    }

    /*
     * Make sure the constructor corrects an invalid direction vector.
     */
    TEST(GeometricObjectsTest, BadNormalDirection)
    {
        Ray ray = Ray(TVector3(10, 120, -11), TVector3(0, 0, 0), 0.04);
        ASSERT_TRUE(VectorsEqual(ray.Direction(), TVector3(0, 0, 1), 0.001));
    }
}
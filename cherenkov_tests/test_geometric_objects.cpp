//
// Created by Matthew Dutson on 1/17/17.
//

#include <gtest/gtest.h>
#include <TMath.h>

#include "geometric_objects.h"
#include "test_helper.h"

using namespace cherenkov_library;
using namespace TMath;

namespace cherenkov_tests
{
    /*
     * Test the default constructor for a plane object.
     */
    TEST(plane, default_construct)
    {
        Plane plane = Plane();
        ASSERT_EQ(plane.Normal(), TVector3(0, 0, 1));
        ASSERT_EQ(plane.Coefficient(), 0);
    }

    /*
     * Check that the normal vector is unit.
     */
    TEST(plane, unit_normal)
    {
        Plane plane = Plane(TVector3(0, 0, 12), TVector3(0, 0, 0));
        ASSERT_EQ(plane.Normal(), TVector3(0, 0, 1));
    }

    /*
     * Test that the correct value of the coefficient is returned.
     */
    TEST(plane, coefficient)
    {
        Plane plane = Plane(TVector3(0, 0, 1), TVector3(45, -12, 33));
        ASSERT_EQ(plane.Coefficient(), 33);
    }

    /*
     * Set a non-default normal vector through the constructor.
     */
    TEST(plane, set_normal)
    {
        Plane plane = Plane(TVector3(1, 1, 1), TVector3(45, -12, 33));
        ASSERT_EQ(plane.Normal(), TVector3(1 / Sqrt(3), 1 / Sqrt(3), 1 / Sqrt(3)));
    }

    /*
     * Make sure the constructor corrects an invalid direction vector.
     */
    TEST(ray, bad_direction)
    {
        Ray ray = Ray(TVector3(10, 120, -11), TVector3(0, 0, 0), 0.04);
        ASSERT_TRUE(VectorsEqual(ray.Direction(), TVector3(0, 0, 1), 0.001));
    }
}
// GeometricTest.cpp
//
// Author: Matthew Dutson
//
// Tests of Geometric.h

#include <gtest/gtest.h>
#include <TMath.h>

#include "Geometric.h"
#include "Helper.h"

using namespace std;
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
    Shower::Params test_params;
    Shower test_shower;

    virtual void SetUp()
    {
        test_ray_1 = Ray(TVector3(7, 2, 3), TVector3(8, 8, 8), 1.7);
        test_ray_2 = Ray(TVector3(0, 0, 0), TVector3(-1, 2, -1), -0.8);

        test_params = Shower::Params();
        test_params.energy = 2.7e19;
        test_params.x_max = 800;
        test_params.n_max = 2.1e10;
        test_params.rho_0 = 0.0012;
        test_params.scale_height = 841300;
        test_params.delta_0 = 0.00029;

        test_shower = Shower(test_params, TVector3(0, 0, 2000000), TVector3(1, -1, -3));
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

    Shower::Params CopyParams()
    {
        return test_params;
    }

    Shower CopyShower()
    {
        return test_shower;
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
     * An invalid_argument exception should be thrown if the plane normal vector is zero.
     */
    TEST_F(GeometricTest, BadNormalDirection)
    {
        try
        {
            Plane(TVector3(0, 0, 0), TVector3(10, 120, -11));
            FAIL() << "Exception not thrown";
        }
        catch(invalid_argument& err)
        {
            ASSERT_EQ(string("Plane normal vector must be nonzero"), err.what());
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
        ASSERT_TRUE(Helper::VectorsEqual(TVector3(1, 1, 1).Unit(), ray.Direction(), 1e-6));
        ASSERT_EQ(1.7, ray.Time());
    }

    /*
     * An invalid_argument exception should be thrown if the Ray direction vector is zero.
     */
    TEST_F(GeometricTest, InvalidDirection)
    {
        try
        {
            Ray(TVector3(1, 2, 3), TVector3(0, 0, 0), 0);
            FAIL() << "Exception not thrown";
        }
        catch(invalid_argument& err)
        {
            ASSERT_EQ(string("Ray direction must be nonzero"), err.what());
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
        ASSERT_TRUE(Helper::VectorsEqual(TVector3(1, 1, 1).Unit(), ray.Direction(), 1e-6));

        ray.Reflect(TVector3(-1, 0, 0));
        ASSERT_TRUE(Helper::VectorsEqual(TVector3(-1, 1, 1).Unit(), ray.Direction(), 1e-6));

        ray.SetDirection(TVector3(-1, 2, 0));
        ASSERT_TRUE(Helper::VectorsEqual(TVector3(-1, 2, 0).Unit(), ray.Direction(), 1e-6));
    }

    /*
     * An invalid_argument exception should be thrown if the Ray direction vector is zero.
     */
    TEST_F(GeometricTest, SetInvalidDirection)
    {
        Ray ray = CopyRay1();
        try
        {
            ray.SetDirection(TVector3(0, 0, 0));
            FAIL() << "Exception not thrown";
        }
        catch(invalid_argument& err)
        {
            ASSERT_EQ(string("Ray direction must be nonzero"), err.what());
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
        ASSERT_TRUE(Helper::VectorsEqual(TVector3(-8, 97, 4), ray.Position(), 1e-6));
        ASSERT_TRUE(Helper::VectorsEqual(TVector3(-15, 95, 1).Unit(), ray.Direction(), 1e-6));
        ASSERT_EQ(1.7 + TVector3(-15, 95, 1).Mag() / Utility::c_cent, ray.Time());
    }

    /*
     * Tests the PropagateToPlane method when the plane is in front of or behind the Ray.
     */
    TEST_F(GeometricTest, PropagateToPlane)
    {
        Ray ray1 = CopyRay1();
        Plane plane1 = Plane(TVector3(1, 0, 0), TVector3(10, 0, 0));
        double init1 = ray1.Time();
        ray1.PropagateToPlane(plane1);
        ASSERT_EQ(ray1.Position().Dot(plane1.Normal()), plane1.Coefficient());
        ASSERT_EQ(TVector3(10, 5, 6), ray1.Position());
        ASSERT_TRUE(ray1.Time() > init1);

        Ray ray2 = CopyRay2();
        Plane plane2 = Plane(TVector3(0, 0, 1), TVector3(0, 0, 2));
        double init2 = ray2.Time();
        ray2.PropagateToPlane(plane2);
        ASSERT_EQ(ray2.Position().Dot(plane2.Normal()), plane2.Coefficient());
        ASSERT_EQ(TVector3(2, -4, 2), ray2.Position());
        ASSERT_TRUE(ray2.Time() < init2);
    }

    /*
     * Tests the PropagateToPlane method when the plane is parallel to the Ray.
     */
    TEST_F(GeometricTest, PropagateToPlaneParallel)
    {
        Ray ray = CopyRay1();
        Plane plane = Plane(TVector3(-2, 1, 1), TVector3(0, 0, 0));
        TVector3 pos_init = ray.Position();
        TVector3 dir_init = ray.Direction();
        double time_init = ray.Time();
        ray.PropagateToPlane(plane);
        ASSERT_TRUE(Helper::VectorsEqual(pos_init, ray.Position(), 1e-6));
        ASSERT_EQ(dir_init, ray.Direction());
        ASSERT_EQ(time_init, ray.Time());
    }

    /*
     * Tests the PlaneImpact method when the plane is in front of or behind the Ray.
     */
    TEST_F(GeometricTest, PlaneImpact)
    {
        Ray ray1 = CopyRay1();
        Plane plane1 = Plane(TVector3(1, 0, 0), TVector3(10, 0, 0));
        ASSERT_EQ(ray1.PlaneImpact(plane1).Dot(plane1.Normal()), plane1.Coefficient());
        ASSERT_EQ(TVector3(10, 5, 6), ray1.PlaneImpact(plane1));

        Ray ray2 = CopyRay2();
        Plane plane2 = Plane(TVector3(0, 0, 1), TVector3(0, 0, 2));
        ASSERT_EQ(ray2.PlaneImpact(plane2).Dot(plane2.Normal()), plane2.Coefficient());
        ASSERT_EQ(TVector3(2, -4, 2), ray2.PlaneImpact(plane2));
    }

    /*
     * Tests the PlaneImpact method when the plane is parallel to the Ray (the Ray's current position should be
     * returned).
     */
    TEST_F(GeometricTest, PlaneImpactParallel)
    {
        Ray ray = CopyRay1();
        Plane plane = Plane(TVector3(-2, 1, 1), TVector3(0, 0, 0));
        ASSERT_EQ(ray.Position(), ray.PlaneImpact(plane));
    }

    /*
     * Test the TimeToPlane method when the plane is in front of or behind the Ray.
     */
    TEST_F(GeometricTest, TimeToPlane)
    {
        Ray ray1 = CopyRay1();
        Plane plane1 = Plane(TVector3(1, 0, 0), TVector3(10, 0, 0));
        ASSERT_EQ(TVector3(3, 3, 3).Mag() / Utility::c_cent, ray1.TimeToPlane(plane1));

        Ray ray2 = CopyRay2();
        Plane plane2 = Plane(TVector3(0, 0, 1), TVector3(0, 0, 2));
        ASSERT_EQ(- TVector3(2, -4, 2).Mag() / Utility::c_cent, ray2.TimeToPlane(plane2));
    }

    /*
     * Tests the TimeToPlane method when the plane is parallel to the Ray (infinity should be returned).
     */
    TEST_F(GeometricTest, TimeToPlaneParallel)
    {
        Ray ray = CopyRay1();
        Plane plane = Plane(TVector3(-2, 1, 1), TVector3(0, 0, 0));
        ASSERT_EQ(Infinity(), ray.TimeToPlane(plane));
    }

    /*
     * Tests the Reflect method. The sign of the normal vector SHOULD NOT matter.
     */
    TEST_F(GeometricTest, Reflect)
    {
        Ray ray1 = CopyRay1();
        ray1.Reflect(TVector3(-1, -1, 0));
        ASSERT_TRUE(Helper::VectorsEqual(TVector3(-1, -1, 1).Unit() , ray1.Direction(), 1e-6));

        Ray ray2 = CopyRay2();
        ray2.Reflect(TVector3(0, 1, 0));
        ASSERT_TRUE(Helper::VectorsEqual(TVector3(-1, -2, -1).Unit() , ray2.Direction(), 1e-6));
    }

    /*
     * An invalid_argument exception should be thrown if normal vector passed to Reflect() is zero.
     */
    TEST_F(GeometricTest, ReflectZeroNormal)
    {
        try
        {
            CopyRay1().Reflect(TVector3(0, 0, 0));
            FAIL() << "Exception not thrown";
        }
        catch(invalid_argument& err)
        {
            ASSERT_EQ(string("Reflection normal vector must be nonzero"), err.what());
        }
    }

    /*
     * Tests a normal application of the Refract method, both from less dense to more dense and vice versa.
     */
    TEST_F(GeometricTest, Refract)
    {
        Ray ray1 = Ray(TVector3(0, 0, 0), TVector3(1, -1, 0), 0.0);
        ray1.Refract(TVector3(0, 1, 0), 1.0, 1.2);
        double theta_1 = ASin(Sin(Pi() / 4.0) / 1.2);
        ASSERT_TRUE(Helper::VectorsEqual(TVector3(Tan(theta_1), -1, 0).Unit(), ray1.Direction(), 1e-6));

        Ray ray2 = Ray(TVector3(0, 0, 0), TVector3(0, -1, 0), 0.0);
        ray2.Refract(TVector3(0, 1, 0), 1.0, 1.7);
        ASSERT_EQ(TVector3(0, -1, 0), ray2.Direction());

        Ray ray3 = Ray(TVector3(0, 0, 0), TVector3(-1, 2, 0), 0.0);
        ray3.Refract(TVector3(0, -1, 0), 1.3, 1.0);
        double theta_3 = ASin(1.3 * Sin(ATan(0.5)));
        ASSERT_TRUE(Helper::VectorsEqual(TVector3(-Tan(theta_3), 1, 0).Unit(), ray3.Direction(), 1e-6));
    }

    /*
     * An invalid_argument exception should be thrown if normal vector passed to Refract() is zero.
     */
    TEST_F(GeometricTest, RefractZeroNormal)
    {
        try
        {
            CopyRay1().Refract(TVector3(0, 0, 0), 1.0, 1.2);
            FAIL() << "Exception not thrown";
        }
        catch(invalid_argument& err)
        {
            ASSERT_EQ(string("Refraction normal vector must be nonzero"), err.what());
        }
    }

    /*
     * Tests the Refract method when total internal reflection occurs (will only be the case if n_out < n_in).
     */
    TEST_F(GeometricTest, RefractCritAngle)
    {
        double n_in = 1.4;
        double theta_c = ASin(1.0 / n_in);

        Ray ray1 = Ray(TVector3(0, 0, 0), TVector3(Tan(theta_c + 0.05), -1, 0), 0.0);
        TVector3 init1 = ray1.Direction();
        ASSERT_FALSE(ray1.Refract(TVector3(0, 1, 0), n_in, 1.0));
        ASSERT_TRUE(Helper::VectorsEqual(init1, ray1.Direction(), 1e-6));

        Ray ray2 = Ray(TVector3(0, 0, 0), TVector3(Tan(theta_c - 0.05), -1, 0), 0.0);
        TVector3 init2 = ray2.Direction();
        ASSERT_TRUE(ray2.Refract(TVector3(0, 1, 0), n_in, 1.0));
        ASSERT_FALSE(Helper::VectorsEqual(init2, ray2.Direction(), 1e-6));
    }

    /*
     * Test the Refract method when an n < 1 is passed. In this case, invalid_argument should be thrown.
     */
    TEST_F(GeometricTest, RefractInvalidN)
    {
        try
        {
            CopyRay1().Refract(TVector3(-1, 0, 0), 0.9, 1.2);
            FAIL() << "Exception not thrown";
        }
        catch(invalid_argument& err)
        {
            ASSERT_EQ(string("Indices of refraction must be at least one"), err.what());
        }
        try
        {
            CopyRay1().Refract(TVector3(0, -1, 0), 1.01, -0.8);
            FAIL() << "Exception not thrown";
        }
        catch(invalid_argument& err)
        {
            ASSERT_EQ(string("Indices of refraction must be at least one"), err.what());
        }
    }

    /*
     * Test the Transform method.
     */
    TEST_F(GeometricTest, Transform)
    {
        Ray ray1 = CopyRay1();
        TRotation rotation1 = TRotation();
        rotation1.RotateZ(Pi() / 6.0);
        ray1.Transform(rotation1);
        double pos_ang = ATan(2.0 / 7.0) + Pi() / 6.0;
        ASSERT_TRUE(Helper::VectorsEqual(TVector3(Sqrt(53.0) * Cos(pos_ang), Sqrt(53.0) * Sin(pos_ang), 3),
                                         ray1.Position(), 1e-6));
        ASSERT_TRUE(Helper::VectorsEqual(TVector3(Sqrt(2.0) * Sin(Pi() / 12.0), Sqrt(2.0) * Cos(Pi() / 12.0), 1).Unit(),
                                         ray1.Direction(), 1e-6));

        Ray ray2 = CopyRay2();
        TRotation rotation2 = TRotation();
        ray2.Transform(rotation2);
        ASSERT_EQ(TVector3(0, 0, 0), ray2.Position());
        ASSERT_TRUE(Helper::VectorsEqual(TVector3(-1, 2, -1).Unit(), ray2.Direction(), 1e-6));
    }

    /*
     * Test the non-default constructor for Shower.
     */
    TEST_F(GeometricTest, UserShower)
    {
        Shower::Params params = CopyParams();
        Shower shower = CopyShower();
        ASSERT_EQ(params.energy, shower.EnergyeV());
        ASSERT_EQ(TVector3(0, 0, 2000000), shower.Position());
        ASSERT_TRUE(Helper::VectorsEqual(TVector3(1, -1, -3).Unit(), shower.Direction(), 1e-6));
    }

    /*
     * An invalid_argument exception should be thrown if a non-positive energy is passed to the Shower constructor.
     */
    TEST_F(GeometricTest, NonPositiveEnergy)
    {
        Shower::Params params = CopyParams();
        params.energy = -1.2e18;
        try
        {
            Shower(params, TVector3(), TVector3(1, 0, 0));
            FAIL() << "Exception not thrown";
        }
        catch(invalid_argument& err)
        {
            ASSERT_EQ(string("Shower energy must be positive"), err.what());
        }
    }

    /*
     * An invalid_argument exception should be thrown if a non-positive x_max is passed to the Shower constructor.
     */
    TEST_F(GeometricTest, NonPositiveXMax)
    {
        Shower::Params params = CopyParams();
        params.x_max = 0.0;
        try
        {
            Shower(params, TVector3(), TVector3(1, 0, 0));
            FAIL() << "Exception not thrown";
        }
        catch(invalid_argument& err)
        {
            ASSERT_EQ(string("Shower XMax must be positive"), err.what());
        }
    }

    /*
     * An invalid_argument exception should be thrown if a non-positive n_max is passed to the Shower constructor.
     */
    TEST_F(GeometricTest, NonPositiveNMax)
    {
        Shower::Params params = CopyParams();
        params.n_max = -0.3;
        try
        {
            Shower(params, TVector3(), TVector3(1, 0, 0));
            FAIL() << "Exception not thrown";
        }
        catch(invalid_argument& err)
        {
            ASSERT_EQ(string("Shower NMax must be positive"), err.what());
        }
    }

    /*
     * An invalid_argument exception should be thrown if a non-positive rho_0 is passed to the Shower constructor.
     */
    TEST_F(GeometricTest, NonPositiveRho0)
    {
        Shower::Params params = CopyParams();
        params.rho_0 = -0.0012;
        try
        {
            Shower(params, TVector3(), TVector3(1, 0, 0));
            FAIL() << "Exception not thrown";
        }
        catch(invalid_argument& err)
        {
            ASSERT_EQ(string("Atmospheric density must be positive"), err.what());
        }
    }

    /*
     * An invalid_argument exception should be thrown if a non-positive scale height is passed to the Shower
     * constructor.
     */
    TEST_F(GeometricTest, NonPositiveScaleH)
    {
        Shower::Params params = CopyParams();
        params.scale_height = -2.3e-19;
        try
        {
            Shower(params, TVector3(), TVector3(1, 0, 0));
            FAIL() << "Exception not thrown";
        }
        catch(invalid_argument& err)
        {
            ASSERT_EQ(string("Scale height must be positive"), err.what());
        }
    }

    /*
     * An invalid_argument exception should be thrown if a non-positive delta_0 is passed to the Shower constructor.
     */
    TEST_F(GeometricTest, NonPositiveDelta0)
    {
        Shower::Params params = CopyParams();
        params.delta_0 = 0.0;
        try
        {
            Shower(params, TVector3(), TVector3(1, 0, 0));
            FAIL() << "Exception not thrown";
        }
        catch(invalid_argument& err)
        {
            ASSERT_EQ(string("Atmospheric delta0 must be positive"), err.what());
        }
    }

    /*
     * Test the Age function. The depth was determined by independently integrating from shower.Position().Z() to
     * infinity.
     */
    TEST_F(GeometricTest, ShowerAge)
    {
        Shower::Params params = CopyParams();
        Shower shower = CopyShower();
        double depth = 93.6905;
        double x = depth / Abs(shower.Direction().CosTheta());
        ASSERT_TRUE(Helper::ValuesEqual(3.0 * x / (x + 2.0 * params.x_max), shower.Age(), 1e-6));
    }

    /*
     * Test the EnergyMeV function.
     */
    TEST_F(GeometricTest, EnergyMeV)
    {
        Shower::Params params = CopyParams();
        Shower shower = CopyShower();
        ASSERT_EQ(params.energy / 10e6, shower.EnergyMeV());
    }

    /*
     * Test the EnergyeV function.
     */
    TEST_F(GeometricTest, EnergyeV)
    {
        Shower::Params params = CopyParams();
        Shower shower = CopyShower();
        ASSERT_EQ(params.energy, shower.EnergyeV());
    }

    /*
     * Test the ImpactParam function. Used http://onlinemschool.com/math/assistance/cartesian_coordinate/p_line/ to
     * calculate this.
     */
    TEST_F(GeometricTest, ImpactParam)
    {
        Shower shower = CopyShower();
        ASSERT_TRUE(Helper::ValuesEqual(852802.87, shower.ImpactParam(), 1e-6));
    }

    /*
     * Test the ImpactAngle function. Impact angle is the angle between the world y-axis and the particle velocity.
     */
    TEST_F(GeometricTest, ImpactAngle)
    {
        Shower shower = CopyShower();
        Plane horiz_plane = Plane(TVector3(0, 0, 1), TVector3());
        ASSERT_EQ(shower.Direction().Angle(shower.PlaneImpact(horiz_plane)), shower.ImpactAngle());
    }

    /*
     * Test the LocalRho function.
     */
    TEST_F(GeometricTest, LocalRho)
    {
        Shower::Params params = CopyParams();
        Shower shower = CopyShower();
        double rho = params.rho_0 * Exp(-shower.Position().Z() / params.scale_height);
        ASSERT_TRUE(Helper::ValuesEqual(rho, shower.LocalRho(), 1e-6));
    }

    /*
     * Test the LocalDelta function.
     */
    TEST_F(GeometricTest, LocalDelta)
    {
        Shower::Params params = CopyParams();
        Shower shower = CopyShower();
        double delta = params.delta_0 * Exp(-shower.Position().Z() / params.scale_height);
        ASSERT_TRUE(Helper::ValuesEqual(delta, shower.LocalDelta(), 1e-6));
    }

    /*
     * Test the GaisserHillas function. The depth was determined by independently integrating from shower.Position().Z()
     * to infinity.
     */
    TEST_F(GeometricTest, GaisserHillas)
    {
        Shower::Params params = CopyParams();
        Shower shower = CopyShower();
        double depth = 93.6905;
        double x = depth / Abs(shower.Direction().CosTheta());
        double term_1 = Power((x + 70.0)/(params.x_max + 70.0), (params.x_max + 70.0) / 70.0);
        double term_2 = Exp((params.x_max - x) / 70.0);
        ASSERT_TRUE(Helper::ValuesEqual(params.n_max * term_1 * term_2, shower.GaisserHillas(), 1e-4));
    }

    /*
     * Test the EThresh function.
     */
    TEST_F(GeometricTest, EThresh)
    {
        Shower shower = CopyShower();
        ASSERT_EQ(Utility::mass_e / Sqrt(2.0 * shower.LocalDelta()), shower.EThresh());
    }

    /*
     * Test the IncrementDepth function. Note that this is defined for zero and negative depths.The depth was determined
     * by independently integrating from shower.Position().Z() to infinity.
     */
    TEST_F(GeometricTest, IncrementDepth)
    {
        Shower::Params params = CopyParams();
        Shower shower = CopyShower();
        double depth = 93.6905;
        double x = depth / Abs(shower.Direction().CosTheta()) + 1.8;
        shower.IncrementDepth(1.8);
        ASSERT_TRUE(Helper::ValuesEqual(3.0 * x / (x + 2.0 * params.x_max), shower.Age(), 1e-3));
    }
}
// Geometric.cpp
//
// Author: Matthew Dutson
//
// Implementation of Geometric.h

#include <TMath.h>

#include "Geometric.h"
#include "Utility.h"

using namespace TMath;

namespace cherenkov_simulator
{
    Plane::Plane() : Plane(TVector3(), TVector3())
    {}

    Plane::Plane(TVector3 normal_vector, TVector3 point)
    {
        // If a zero normal vector is passed, use (0, 0, 1) instead.
        normal = (normal_vector == TVector3(0, 0, 0)) ? TVector3(0, 0, 1) : normal_vector.Unit();

        // d = a * x_0 + b * y_0 + c * z_0
        coefficient = normal.Dot(point);
    }

    TVector3 Plane::Normal()
    {
        return normal;
    }

    double Plane::Coefficient()
    {
        return coefficient;
    }

    bool Plane::InFrontOf(TVector3 direction)
    {
        Ray outward_ray = Ray(TVector3(), direction, 0);
        return outward_ray.TimeToPlane(*this) > 0;
    }

    Ray::Ray(TVector3 position, TVector3 direction, double time)
    {
        // If a zero direction vector is passed, use (0, 0, 1) instead.
        SetDirection((direction == TVector3(0, 0, 0)) ? TVector3(0, 0, 1) : direction);

        // All times and starting positions are assumed valid.
        this->time = time;
        this->position = position;
    }

    TVector3 Ray::Position()
    {
        return position;
    }

    TVector3 Ray::Velocity()
    {
        return velocity;
    }

    TVector3 Ray::Direction()
    {
        return velocity.Unit();
    }

    void Ray::SetDirection(TVector3 direction)
    {
        // Use the speed of light in centimeters/second.
        velocity = direction.Unit() * Utility::c_cent;
    }

    double Ray::Time()
    {
        return time;
    }

    void Ray::PropagateToPoint(TVector3 destination)
    {
        // If needed, change the direction so the ray is pointing toward the destination.
        TVector3 displacement = destination - position;
        SetDirection(displacement);

        // Move the ray forward until it reaches the destination.
        IncrementPosition(displacement.Mag());
    }

    void Ray::PropagateToPlane(Plane plane)
    {
        PropagateToPoint(PlaneImpact(plane));
    }

    TVector3 Ray::PlaneImpact(Plane plane)
    {
        // If the ray and the plane are exactly parallel, return the current position of the ray.
        double time = TimeToPlane(plane);
        return (time == Infinity()) ? position : position + time * velocity;
    }

    double Ray::TimeToPlane(Plane plane)
    {
        TVector3 normal = plane.Normal();
        double coefficient = plane.Coefficient();

        // Check the edge case where the vector is perfectly parallel to the plane.
        if (normal.Dot(velocity) == 0) return Infinity();
        return (coefficient - normal.Dot(position)) / normal.Dot(velocity);
    }

    void Ray::Reflect(TVector3 normal)
    {
        velocity -= 2 * velocity.Dot(normal.Unit()) * normal.Unit();
    }

    bool Ray::Refract(TVector3 normal, double n_in, double n_out)
    {
        // Reverse the normal vector so it points in the direction of the incoming ray
        normal = -normal;

        double angle_in = velocity.Angle(normal);

        // If the current velocity and normal vector are parallel, don't do anything.
        if (angle_in == 0) return true;

        // If we're more than 90 degrees from the normal, we're coming from the wrong side of the lens (reversed normal
        // should point in the same direction as the incoming ray)
        if (angle_in > PiOver2()) return false;

        // Refract and rotate around some vector perpendicular to both the ray and plane normal.
        double angle_out = ASin(n_in * Sin(angle_in) / n_out);
        TVector3 mutual_norm = normal.Cross(velocity);
        velocity.Rotate(angle_out - angle_in, mutual_norm);
        return true;
    }

    void Ray::Transform(TRotation rotation)
    {
        velocity = rotation * velocity;
        position = rotation * position;
    }

    void Ray::IncrementPosition(double distance)
    {
        double time = distance / Utility::c_cent;
        IncrementTime(time);
    }

    void Ray::IncrementTime(double time_step)
    {
        time += time_step;
        position += time_step * velocity;
    }

    Shower::Shower(Params params, TVector3 position, TVector3 direction, double time) : Ray(position, direction, time)
    {
        this->energy = params.energy;
        this->x_max = params.x_max;
        this->n_max = params.n_max;
        this->rho_0 = params.rho_0;
        this->scale_height = params.scale_height;
        this->delta_0 = params.delta_0;
    }

    double Shower::Age()
    {
        return 3.0 * X() / (X() + 2.0 * x_max);
    }

    double Shower::EnergyMeV()
    {
        return energy / (10.0e6);
    }

    double Shower::LocalRho()
    {
        return rho_0 * Exp(-position.Z() / scale_height);
    }

    double Shower::LocalDelta()
    {
        return delta_0 * Exp(-position.Z() / scale_height);
    }

    double Shower::GaisserHillas()
    {
        double pow = Power((X() - x_0) / (x_max - x_0), (x_max - x_0) / gh_lambda);
        double exp = Exp((x_max - X()) / gh_lambda);
        return n_max * pow * exp;
    }

    double Shower::EThresh()
    {
        return Utility::mass_e / Sqrt(2 * LocalDelta());
    }

    void Shower::IncrementDepth(double depth)
    {
        // We make the simplifying assumption that the atmospheric density is constant over a single step
        IncrementPosition(depth / LocalRho());
    }

    double Shower::X()
    {
        // See 1/25 depth integration notes.
        double cos_theta = Abs(velocity.CosTheta());
        return scale_height * rho_0 / cos_theta * Exp(-position.Z() / scale_height);
    }
}

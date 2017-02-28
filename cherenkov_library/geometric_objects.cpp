// geometric_objects.cpp
// cherenkov_library
//
// Created by Matthew Dutson on 9/8/16.
//
// Contains the implementation of methods in geometric_objects.h.

#include "geometric_objects.h"
#include "utility.h"

using namespace TMath;

namespace cherenkov_library
{
    Plane::Plane() : Plane(TVector3(), TVector3())
    {}

    Plane::Plane(TVector3 normal_vector, TVector3 point)
    {
        // If a zero normal vector is passed, use (0, 0, 1) instead.
        if (normal_vector == TVector3(0, 0, 0))
        {
            normal = TVector3(0, 0, 1);
        }
        else
        {
            normal = normal_vector.Unit();
        }

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
        if (direction == TVector3(0, 0, 0))
        {
            SetDirection(TVector3(0, 0, 1));
        }
        else
        {
            SetDirection(direction);
        }

        // All times and starting positions are assumed valid.
        current_time = time;
        current_position = position;
    }

    double Ray::Time()
    {
        return current_time;
    }

    TVector3 Ray::Position()
    {
        return current_position;
    }

    TVector3 Ray::Velocity()
    {
        return current_velocity;
    }

    TVector3 Ray::Direction()
    {
        return current_velocity.Unit();
    }

    void Ray::Transform(TRotation rotation)
    {
        current_velocity = rotation * current_velocity;
        current_position = rotation * current_position;
    }

    void Ray::SetDirection(TVector3 direction)
    {
        // Use the speed of light in centimeters/second.
        current_velocity = direction.Unit() * CentC();
    }

    void Ray::PropagateToPoint(TVector3 destination)
    {
        // If needed, change the direction so the ray is pointing toward the destination.
        TVector3 displacement = destination - current_position;
        SetDirection(displacement);

        // Move the ray forward until it reaches the destination.
        IncrementPosition(displacement.Mag());
    }

    void Ray::PropagateToPlane(Plane plane)
    {
        PropagateToPoint(PlaneImpact(plane));
    }

    double Ray::TimeToPlane(Plane plane)
    {
        TVector3 normal = plane.Normal();
        double coefficient = plane.Coefficient();

        // Check the edge case where the vector is perfectly parallel to the plane.
        if (normal.Dot(current_velocity) == 0)
        {
            return Infinity();
        }
        else
        {
            return (coefficient - normal.Dot(current_position)) / normal.Dot(current_velocity);
        }
    }

    TVector3 Ray::PlaneImpact(Plane plane)
    {
        double time = TimeToPlane(plane);

        // If the ray and the plane are exactly parallel, return the current position of the ray.
        if (time == Infinity())
        {
            return current_position;
        }
        else
        {
            return current_position + time * current_velocity;
        }
    }

    void Ray::Reflect(TVector3 normal)
    {
        current_velocity -= 2 * current_velocity.Dot(normal.Unit()) * normal.Unit();
    }

    bool Ray::Refract(TVector3 normal, double n_in, double n_out)
    {
        double angle_in = current_velocity.Angle(normal);

        // If the current velocity and normal vector are parallel, don't do anything.
        if (angle_in == 0)
        {
            return true;
        }

            // If we're more than 90 degrees from the normal, we're coming from the wrong side of the lens. This can
            // occur when the shower is nearly 90 degrees off the detector axis.
        else if (angle_in > PiOver2())
        {
            return false;
        }

            // Refract and rotate around some vector perpendicular to both the ray and plane normal.
        else
        {
            double angle_out = ASin(n_in * Sin(angle_in) / n_out);
            TVector3 mutual_norm = normal.Cross(current_velocity);
            current_velocity.Rotate(angle_out - angle_in, mutual_norm);
            return true;
        }
    }

    void Ray::IncrementPosition(double distance)
    {
        double time = distance / CentC();
        IncrementTime(time);
    }

    void Ray::IncrementTime(double time)
    {
        current_time += time;
        current_position += time * current_velocity;
    }

    Shower::Shower(Params params, TVector3 position, TVector3 direction, double time) : Ray(position, direction, time)
    {
        start_position = position;
        this->energy = params.energy;
        this->x_0 = params.x_0;
        this->x_max = params.x_max;
        this->n_max = params.n_max;
        this->rho_0 = params.rho_0;
        this->scale_height = params.scale_height;
        this->delta_0 = params.delta_0;
    }

    double Shower::Age()
    {
        return 3.0 * X() / (X() + 2.0 * XMax());
    }

    double Shower::EnergyMeV()
    {
        return energy / (10.0e6);
    }

    double Shower::X()
    {
        // See 1/25 depth integration notes.
        double cos_theta = Abs(current_position.CosTheta());
        return scale_height * rho_0 / cos_theta * Exp(-current_position.Z() / scale_height);
    }

    double Shower::X0()
    {
        return x_0;
    }

    double Shower::XMax()
    {
        return x_max;
    }

    double Shower::NMax()
    {
        return n_max;
    }

    double Shower::LocalRho()
    {
        return rho_0 * Exp(-current_position.Z() / scale_height);
    }

    double Shower::LocalDelta()
    {
        return delta_0 * Exp(-current_position.Z() / scale_height);
    }

    void Shower::IncrementDepth(double depth)
    {
        // We make the simplifying assumption that the density of the atmosphere is approximately constant in the
        // range of a single step.
        // TODO: Use an exact formula (the integral of the atmospheric density is easy to find).
        IncrementPosition(depth / LocalRho());
    }
}

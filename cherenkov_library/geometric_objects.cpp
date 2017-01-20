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
        current_velocity -= 2 * current_velocity.Dot(normal) * normal;
    }

    void Ray::Refract(TVector3 normal, double n_in, double n_out)
    {
        // If the current velocity and normal vector are parallel, don't do anything.
        if (normal.Cross(current_velocity) != TVector3(0, 0, 0))
        {
            // Apply Snell's law.
            double angle_in = current_velocity.Angle(normal);
            double angle_out = ASin(n_in * Sin(angle_in) / n_out);

            // Refract and rotate around some vector perpendicular to both the ray and plane normal.
            TVector3 mutual_norm = normal.Cross(current_velocity);
            current_velocity.Rotate(angle_out - angle_in, mutual_norm);
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

    double Shower::X()
    {
        // TODO: Should we be integrating from the shower's starting point or from the top of the atmosphere?
        TVector3 point1 = Position();
        TVector3 point2 = start_position;
        TVector3 diff = point2 - point1;
        double cosine = Abs(diff.Z() / diff.Mag());
        double vertical_depth =
                -rho_0 / scale_height * (Exp(-point2.Z() / scale_height) - Exp(-point1.Z() / scale_height));
        return vertical_depth / cosine;
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
        double vertical_distance =
                -scale_height * Log(Exp(-current_position.Z() / scale_height) + scale_height * depth / rho_0);
        double total_distance = vertical_distance / Abs(current_velocity.CosTheta());
        IncrementPosition(total_distance);
    }
}

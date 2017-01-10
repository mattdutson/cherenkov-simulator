// geometric_objects.cpp
// cherenkov_library
//
// Created by Matthew Dutson on 9/8/16.
//
//

#include "geometric_objects.h"
#include "utility.h"
#include "TMath.h"

using namespace TMath;

namespace cherenkov_simulator
{
    Ray::Ray(double time, TVector3 position, TVector3 direction)
    {
        // TODO: Define behavior when the direction vector is zero
        current_time = time;
        current_position = position;
        SetDirection(direction);
    }

    void Ray::IncrementTime(double time)
    {
        current_time += time;
        current_position += time * current_velocity;
    }

    double Ray::TimeToPlane(Plane plane)
    {
        TVector3 normal = plane.Normal();
        double coefficient = plane.Coefficient();

        // TODO: Determine behavior when normal.Dot(current_velocity) is zero and behavior when time is negative
        return (coefficient - normal.Dot(current_position)) / normal.Dot(current_velocity);
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

    void Ray::SetDirection(TVector3 direction)
    {
        current_velocity = direction.Unit() * C();
    }

    void Ray::Reflect(TVector3 normal)
    {
        current_velocity -= 2 * current_velocity.Dot(normal) * normal;
    }

    void Ray::Refract(TVector3 normal, double n_in, double n_out)
    {
        TVector3 mutual_norm = normal.Cross(current_velocity).Unit();
        double angle_in = current_velocity.Angle(normal);

        // Now refract using Snell's Law.
        double angle_out = ASin(n_in * Sin(angle_in) / n_out);

        current_velocity.Rotate(angle_out - angle_in, mutual_norm);
    }

    void Ray::PropagateToPoint(TVector3 point)
    {
        TVector3 displacement = point - current_position;
        SetDirection(displacement);
        double time = displacement.Mag() / C();
        IncrementTime(time);
    }

    void Ray::PropagateToPlane(Plane plane)
    {
        double time = TimeToPlane(plane);
        IncrementTime(time);
    }

    Shower::Shower(TVector3 position, TVector3 direction, double x_0, double x_max, double n_max, double time) : Ray(
            time, position, direction)
    {
        start_position = position;
        this->x_0 = x_0;
        this->x_max = x_max;
        this->n_max = n_max;
    }

    TVector3 Shower::StartPosition()
    {
        return start_position;
    }

    TVector3 Shower::GroundImpact()
    {
        return ground_impact;
    }

    double Shower::X0()
    {
        return x_0;
    }

    double Shower::X()
    {
        TVector3 point1 = Position();
        TVector3 point2 = StartPosition();
        TVector3 diff = point2 - point1;
        double cosine = Abs(diff.Z() / diff.Mag());
        double vertical_depth =
                -rho_0 / scale_height * (Exp(-point2.Z() / scale_height) - Exp(-point1.Z() / scale_height));
        return vertical_depth / cosine;
    }

    double Shower::XMax()
    {
        return x_max;
    }

    double Shower::NMax()
    {
        return n_max;
    }

    double Shower::Age()
    {
        double x = X();
        return 3.0 * x / (x + 2.0 * XMax());
    }

    Plane::Plane(TVector3 normal_vector, TVector3 point)
    {
        normal = normal_vector.Unit();
        coefficient = normal.Dot(point);
    }

    double Shower::IncrementDepth(double depth)
    {
        double vertical_distance =
                -scale_height * Log(Exp(-current_position.Z() / scale_height) + scale_height * depth / rho_0);
        double total_distance = vertical_distance / Abs(current_velocity.CosTheta());
        IncrementPosition(total_distance);
        return total_distance;
    }
}

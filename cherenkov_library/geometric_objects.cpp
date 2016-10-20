// geometric_objects.cpp
// cherenkov_library
//
// Created by Matthew Dutson on 9/8/16.
//
//

#include <ImageIO/ImageIO.h>
#include "geometric_objects.h"
#include "common.h"

namespace cherenkov_simulator
{

    double ConstantIntensity::GetIntensity(Shower shower)
    {
        return 0;
    }

    Ray::Ray(double time, TVector3 position, TVector3 direction)
    {
        // TODO: Define behavior when the direction vector is zero
        current_time = time;
        current_position = position;
        SetDirection(direction);
    }

    void Ray::IncrementPosition(double time)
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
        current_velocity = direction.Unit() * light_speed;
    }

    void Ray::Reflect(TVector3 normal)
    {
        current_velocity -= 2 * current_velocity.Dot(normal) * normal;
    }

    void Ray::PropagateToPoint(TVector3 point)
    {
        TVector3 displacement = point - current_position;
        SetDirection(displacement);
        double time = displacement.Mag() / light_speed;
        IncrementPosition(time);
    }

    void Ray::PropagateToPlane(Plane plane)
    {
        double time = TimeToPlane(plane);
        IncrementPosition(time);
    }

    Shower::Shower(double time, TVector3 position, TVector3 direction, IntensityFunctor func) : Ray(time, position,
                                                                                                    direction)
    {
        start_time = time;
        start_position = position;
        intensity_functor = func;
    }

    // TODO: Implement this method
    int Shower::NumberFluorescencePhotons()
    {
        return 0;
    }

    // TODO: Implement this method
    int Shower::NumberCherenkovPhotons()
    {
        return 0;
    }

    // TODO: Implement this method
    Ray Shower::GenerateCherenkovPhoton()
    {
        return Ray(0, TVector3(), TVector3());
    }

    Plane::Plane(TVector3 normal, TVector3 point)
    {
        normal_vector = normal.Unit();
        coefficient = normal.Dot(point);
    }
}

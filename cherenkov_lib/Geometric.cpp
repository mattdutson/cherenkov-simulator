// Geometric.cpp
//
// Author: Matthew Dutson
//
// Implementation of Geometric.h

#include <TMath.h>

#include "Geometric.h"

using namespace TMath;
using namespace std;

namespace cherenkov_simulator
{
    Plane::Plane() : Plane(TVector3(0, 0, 1), TVector3())
    {}

    Plane::Plane(TVector3 normal_vector, TVector3 point)
    {
        if (normal_vector.Mag2() == 0) throw std::runtime_error("Plane normal vector must be nonzero");
        normal = normal_vector;
        coefficient = normal.Dot(point);
    }

    TVector3 Plane::Normal() const
    {
        return normal;
    }

    double Plane::Coefficient() const
    {
        return coefficient;
    }

    bool Plane::InFrontOf(TVector3 direction) const
    {
        Ray outward_ray = Ray(TVector3(), direction, 0);
        return outward_ray.TimeToPlane(*this) > 0;
    }

    Ray::Ray() : Ray(TVector3(), TVector3(0, 0, 1), 0)
    {}

    Ray::Ray(TVector3 position, TVector3 direction, double time)
    {
        if (direction.Mag2() == 0) throw std::runtime_error("Ray direction must be nonzero");
        SetDirection(direction);
        this->time = time;
        this->position = position;
    }

    TVector3 Ray::Position() const
    {
        return position;
    }

    TVector3 Ray::Velocity() const
    {
        return velocity;
    }

    TVector3 Ray::Direction() const
    {
        return velocity.Unit();
    }

    void Ray::SetDirection(TVector3 direction)
    {
        // Use the speed of light in centimeters/second.
        velocity = direction.Unit() * Utility::c_cent;
    }

    double Ray::Time() const
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
        PropagateToPoint(PlaneImpact(std::move(plane)));
    }

    TVector3 Ray::PlaneImpact(Plane plane) const
    {
        // If the ray and the plane are exactly parallel, return the current position of the ray.
        double time = TimeToPlane(std::move(plane));
        return (time == Infinity()) ? position : position + time * velocity;
    }

    double Ray::TimeToPlane(Plane plane) const
    {
        TVector3 normal = plane.Normal();
        double coefficient = plane.Coefficient();

        // Check the edge case where the vector is perfectly parallel to the plane.
        if (normal.Dot(velocity) == 0) return Infinity();
        return (coefficient - normal.Dot(position)) / normal.Dot(velocity);
    }

    void Ray::Reflect(TVector3 normal)
    {
        // TODO: Would it be a good idea to update this with the SetDirection method?
        velocity -= 2 * velocity.Dot(normal.Unit()) * normal.Unit();
    }

    bool Ray::Refract(TVector3 normal, double n_in, double n_out)
    {
        // Reverse the normal vector so it points in the direction of the incoming ray
        // TODO: Perform a check and update here on the direction.
        // TODO: Perform a check on the critical angle of the substance.
        // TODO: Would it be a good idea to update the velocity with the SetDirection method?
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

    Shower::Shower() : Ray()
    {}

    Shower::Shower(Params params, TVector3 position, TVector3 direction, double time) : Ray(position, direction, time)
    {
        this->energy = params.energy;
        this->x_max = params.x_max;
        this->n_max = params.n_max;
        this->rho_0 = params.rho_0;
        this->scale_height = params.scale_height;
        this->delta_0 = params.delta_0;
    }

    double Shower::Age() const
    {
        return 3.0 * X() / (X() + 2.0 * x_max);
    }

    double Shower::EnergyMeV() const
    {
        return energy / (10.0e6);
    }

    double Shower::EnergyeV() const
    {
        return energy;
    }

    double Shower::ImpactParam() const
    {
        // See http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
        TVector3 point1 = Position();
        TVector3 point2 = Position() + Direction();
        return (point1.Cross(point2)).Mag();
    }

    double Shower::ImpactAngle() const
    {
        return Direction().Angle(TVector3(0, 1, 0));
    }

    double Shower::LocalRho() const
    {
        return rho_0 * Exp(-position.Z() / scale_height);
    }

    double Shower::LocalDelta() const
    {
        return delta_0 * Exp(-position.Z() / scale_height);
    }

    double Shower::GaisserHillas() const
    {
        double pow = Power((X() - x_0) / (x_max - x_0), (x_max - x_0) / gh_lambda);
        double exp = Exp((x_max - X()) / gh_lambda);
        return n_max * pow * exp;
    }

    double Shower::EThresh() const
    {
        return Utility::mass_e / Sqrt(2 * LocalDelta());
    }

    void Shower::IncrementDepth(double depth)
    {
        // We make the simplifying assumption that the atmospheric density is constant over a single step
        IncrementPosition(depth / LocalRho());
    }

    string Shower::Header()
    {
        return "Psi, Impact";
    }

    string Shower::ToString() const
    {
        return to_string(ImpactAngle()) + ", " + Utility::KmString(ImpactParam());
    }

    double Shower::X() const
    {
        // See 1/25 depth integration notes.
        double cos_theta = Abs(velocity.CosTheta());
        return scale_height * rho_0 / cos_theta * Exp(-position.Z() / scale_height);
    }
}

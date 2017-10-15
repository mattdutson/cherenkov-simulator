// Geometric.cpp
//
// Author: Matthew Dutson
//
// Implementation of Geometric.h

#include <TMath.h>

#include "Geometric.h"

using namespace std;
using namespace TMath;

namespace cherenkov_simulator
{
    Plane::Plane() : Plane(TVector3(0, 0, 1), TVector3())
    {}

    Plane::Plane(TVector3 normal, TVector3 point)
    {
        if (normal.Mag2() == 0) 
            throw invalid_argument("Plane normal vector must be nonzero");
        this->normal = normal.Unit();
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
        if (direction.Mag2() == 0)
            throw invalid_argument("Ray direction must be nonzero");
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
        double time = TimeToPlane(move(plane));
        if (time != Infinity())
            IncrementTime(time);
    }

    TVector3 Ray::PlaneImpact(Plane plane) const
    {
        double time = TimeToPlane(move(plane));
        if (time != Infinity())
            return position + time * velocity;
        return position;
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
        if(normal.Mag2() == 0)
            throw invalid_argument("Reflection normal vector must be nonzero");
        velocity -= 2 * velocity.Dot(normal.Unit()) * normal.Unit();
    }

    bool Ray::Refract(TVector3 normal, double n_in, double n_out)
    {
        // Reverse the normal vector so it points in the direction of the incoming ray
        // TODO: Perform a check and update here on the direction.
        // TODO: Perform a check on the critical angle of the substance.
        // TODO: Would it be a good idea to update the velocity with the SetDirection method?
        if(normal.Mag2() == 0)
            throw invalid_argument("Refraction normal vector must be nonzero");
        if (n_in < 1 || n_out < 1)
            throw invalid_argument("Indices of refraction must be at least one");
        normal = -normal;

        double angle_in = velocity.Angle(normal);

        // If the current velocity and normal vector are parallel, don't do anything.
        if (angle_in == 0) return true;

        // If we're more than 90 degrees from the normal, we're coming from the wrong side of the lens (reversed normal
        // should point in the same direction as the incoming ray)
        double theta_c = ASin(n_out / n_in);
        if (angle_in > PiOver2() || angle_in > theta_c) return false;

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

    Shower::Shower() : Ray() {}

    Shower::Shower(Params params, TVector3 position, TVector3 direction, double time) : Ray(position, direction, time)
    {
        this->energy = params.energy;
        this->x_max = params.x_max;
        this->n_max = params.n_max;
        this->rho_0 = params.rho_0;
        this->scale_height = params.scale_height;
        this->delta_0 = params.delta_0;

        if (energy <= 0.0)
            throw invalid_argument("Shower energy must be positive");
        if (x_max <= 0.0)
            throw invalid_argument("Shower XMax must be positive");
        if (n_max <= 0.0)
            throw invalid_argument("Shower NMax must be positive");
        if (rho_0 <= 0.0)
            throw invalid_argument("Atmospheric density must be positive");
        if (scale_height <= 0.0)
            throw invalid_argument("Scale height must be positive");
        if (delta_0 <= 0.0)
            throw invalid_argument("Atmospheric delta0 must be positive");
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
        Plane horiz_plane = Plane(TVector3(0, 0, 1), TVector3(0, 0, 0));
        return Direction().Angle(PlaneImpact(horiz_plane));
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

// Geometric.cpp
//
// Author: Matthew Dutson
//
// Implementation of Geometric.h.

#include <TMath.h>

#include "Geometric.h"

using namespace std;
using namespace TMath;

namespace cherenkov_simulator
{
    Plane::Plane() : Plane(TVector3(0, 0, 1), TVector3()) {}

    Plane::Plane(TVector3 normal, TVector3 point)
    {
        if (normal.Mag2() == 0) 
            throw invalid_argument("Plane normal vector must be nonzero");
        this->normal = normal.Unit();
        this->coeff = normal.Dot(point);
    }

    TVector3 Plane::Normal() const
    {
        return normal;
    }

    double Plane::Coefficient() const
    {
        return coeff;
    }

    bool Plane::InFrontOf(TVector3 direction) const
    {
        Ray outward_ray = Ray(TVector3(), direction, 0);
        return outward_ray.TimeToPlane(*this) > 0;
    }

    Ray::Ray() : Ray(TVector3(), TVector3(0, 0, 1), 0) {}

    Ray::Ray(TVector3 position, TVector3 direction, double cur_time)
    {
        SetDirection(direction);
        this->cur_time = cur_time;
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
        return cur_time;
    }

    void Ray::PropagateToPoint(TVector3 destination)
    {
        TVector3 displacement = destination - position;
        SetDirection(displacement);
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
        if (normal.Dot(velocity) == 0)
            return Infinity();
        return (plane.Coefficient() - normal.Dot(position)) / normal.Dot(velocity);
    }

    void Ray::Reflect(TVector3 normal)
    {
        if(normal.Mag2() == 0)
            throw invalid_argument("Reflection normal vector must be nonzero");
        SetDirection(velocity - 2 * velocity.Dot(normal.Unit()) * normal.Unit());
    }

    bool Ray::Refract(TVector3 normal, double n_in, double n_out)
    {
        if(normal.Mag2() == 0)
            throw invalid_argument("Refraction normal vector must be nonzero");
        if (n_in < 1 || n_out < 1)
            throw invalid_argument("Indices of refraction must be at least one");

        normal = -normal;
        double angle_in = velocity.Angle(normal);
        if (angle_in > PiOver2() || angle_in > ASin(n_out / n_in)) return false;

        if (angle_in == 0) return true;
        double angle_out = ASin(n_in * Sin(angle_in) / n_out);
        velocity.Rotate(angle_out - angle_in, normal.Cross(velocity));
        return true;
    }

    void Ray::Transform(TRotation rotation)
    {
        velocity *= rotation;
        position *= rotation;
    }

    void Ray::IncrementPosition(double distance)
    {
        IncrementTime(distance / Utility::c_cent);
    }

    void Ray::IncrementTime(double time_step)
    {
        cur_time += time_step;
        position += time_step * velocity;
    }

    Shower::Params::Params()
    {
        energ = 1.0;
        x_max = 1.0;
        n_max = 1.0;
        rho_0 = 1.0;
        atm_h = 1.0;
        del_0 = 1.0;
    }

    Shower::Shower() : Ray() {}

    Shower::Shower(Params params, TVector3 position, TVector3 direction, double time) : Ray(position, direction, time)
    {
        energ = params.energ;
        x_max = params.x_max;
        n_max = params.n_max;
        rho_0 = params.rho_0;
        atm_h = params.atm_h;
        del_0 = params.del_0;

        if (energ <= 0.0)
            throw invalid_argument("Shower energy must be positive");
        if (x_max <= 0.0)
            throw invalid_argument("Shower XMax must be positive");
        if (n_max <= 0.0)
            throw invalid_argument("Shower NMax must be positive");
        if (rho_0 <= 0.0)
            throw invalid_argument("Atmospheric density must be positive");
        if (atm_h <= 0.0)
            throw invalid_argument("Scale height must be positive");
        if (del_0 <= 0.0)
            throw invalid_argument("Atmospheric delta0 must be positive");
    }

    double Shower::Age() const
    {
        return 3.0 * X() / (X() + 2.0 * x_max);
    }

    double Shower::EnergyMeV() const
    {
        return energ / (10.0e6);
    }

    double Shower::EnergyeV() const
    {
        return energ;
    }

    double Shower::ImpactParam() const
    {
        return Position().Cross(Position() + Direction()).Mag();
    }

    double Shower::ImpactAngle() const
    {
        return Direction().Angle(PlaneImpact(Plane()));
    }

    double Shower::LocalRho() const
    {
        return rho_0 * Exp(-position.Z() / atm_h);
    }

    double Shower::LocalDelta() const
    {
        return del_0 * Exp(-position.Z() / atm_h);
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
        return atm_h * rho_0 / Abs(velocity.CosTheta()) * Exp(-position.Z() / atm_h);
    }
}

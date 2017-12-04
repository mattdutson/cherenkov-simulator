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
        velocity = direction.Unit() * c_cent;
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
        IncrementTime(distance / c_cent);
    }

    void Ray::IncrementTime(double time_step)
    {
        cur_time += time_step;
        position += time_step * velocity;
    }

    Shower::Shower() : Ray() {}

    Shower::Shower(double energy, double elevation, TVector3 position, TVector3 direction, double time) : Ray(position, direction, time)
    {
        this->energy = energy;
        this->elevation = elevation;

        if (energy <= 0.0)
            throw invalid_argument("Shower energy must be positive");
    }

    double Shower::Age() const
    {
        return 3.0 * X() / (X() + 2.0 * XMax());
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
        return Position().Cross(Position() + Direction()).Mag();
    }

    double Shower::ImpactAngle() const
    {
        return Direction().Angle(PlaneImpact(Plane()));
    }

    double Shower::LocalRho() const
    {
        return rho_sea * Exp(-(position.Z() + elevation) / scale_h);
    }

    double Shower::LocalDelta() const
    {
        return (ref_sea - 1.0) * Exp(-(position.Z() + elevation) / scale_h);
    }

    double Shower::GaisserHillas() const
    {
        double pow = Power((X() - x_0) / (XMax() - x_0), (XMax() - x_0) / gh_lambda);
        double exp = Exp((XMax() - X()) / gh_lambda);
        return NMax() * pow * exp;
    }

    double Shower::EThresh() const
    {
        return mass_e / Sqrt(2 * LocalDelta());
    }

    void Shower::IncrementDepth(double depth)
    {
        IncrementPosition(depth / LocalRho());
    }

    string Shower::Header()
    {
        return "Angle(deg),Impact(km),Ground(km)";
    }

    string Shower::ToString(Plane ground_plane) const
    {
        return to_string(ImpactAngle() * 180.0 / Pi()) + "," + Utility::KmString(ImpactParam())
               + "," + Utility::KmString(PlaneImpact(move(ground_plane)).Mag());
    }

    double Shower::X() const
    {
        return scale_h * rho_sea / Abs(velocity.CosTheta()) * Exp(-(position.Z() + elevation) / scale_h);
    }

    double Shower::XMax() const
    {
        return x_max_1 + x_max_2 * (Log10(energy) - x_max_3);

    }

    double Shower::NMax() const
    {
        return energy / n_ratio;
    }
}

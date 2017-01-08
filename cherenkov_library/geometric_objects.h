// geometric_objects.h
// cherenkov_library
//
// Created by Matthew Dutson on 9/8/16.
//
//

#ifndef geometric_objects_h
#define geometric_objects_h

#include "TVector3.h"
#include "TRotation.h"

namespace cherenkov_simulator
{
    class Plane
    {
    private:

        TVector3 normal;

        double coefficient;

    public:

        Plane(TVector3 normal, TVector3 point);

        TVector3 Normal();

        double Coefficient();
    };

    class Ray
    {
    protected:

        double current_time;

        TVector3 current_position;

        TVector3 current_velocity;

    public:

        Ray(double time, TVector3 position, TVector3 direction);

        void IncrementTime(double time);

        void IncrementPosition(double position);

        double TimeToPlane(Plane p);

        double Time();

        TVector3 Position();

        TVector3 Velocity();

        void SetDirection(TVector3 direction);

        void Reflect(TVector3 normal);

        /*
         * Refracts the ray across the normal vector using the incident and outward indices of refraction specified.
         */
        void Refract(TVector3 normal, double n_in, double n_out);

        void PropagateToPoint(TVector3 point);

        void PropagateToPlane(Plane plane);

        void Transform(TRotation t);
    };

    class Shower : public Ray
    {
    private:

        TVector3 start_position;

        TVector3 ground_impact;

        double x_0;

        double x_max;

        double n_max;

        // The density of the atmosphere at the origin
        double rho_0;

        // The atmospheric scale height
        double scale_height;

    public:

        Shower(TVector3 position, TVector3 direction, double x_0, double x_max, double n_max, double time = 0);

        void Transform(TRotation t);

        /*
         * Returns the starting position of the shower.
         */
        TVector3 StartPosition();

        /*
         * Returns the impact point of the shower.
         */
        TVector3 GroundImpact();

        /*
         * Returns the slant depth of the first shower interaction.
         */
        double X0();

        /*
         * Returns the slant depth of the shower maximum.
         */
        double XMax();

        /*
         * Returns the number of electrons at the shower maximum.
         */
        double NMax();
    };
}

#endif

// Geometric.h
//
// Author: Matthew Dutson
//
// Contains the definitions of the Plane, Shower, and Ray classes.

#ifndef GEOMETRIC_H
#define GEOMETRIC_H

#include <TRotation.h>
#include <TVector3.h>

namespace cherenkov_simulator
{
    /*
     * A class representing a plane in 3D space.
     */
    class Plane
    {
    public:

        /*
         * The default constructor. Passes two zero vectors to the two-parameter constructor. The behavior for dealing
         * with a zero normal vector is defined in the other constructor.
         */
        Plane();

        /*
         * Constructs a plane with the specified normal vector and with a coefficient "d" which causes the input point
         * to lie on the plane. If (0, 0, 0) is passed as the normal vector, (0, 0, 1) is used instead.
         */
        Plane(TVector3 normal, TVector3 point);

        /*
         * Returns a unit vector normal to the plane.
         */
        TVector3 Normal();

        /*
         * Returns the coefficient "d" in the plane equation.
         */
        double Coefficient();

        /*
         * Returns true if a ray going outward from the origin in the specified direction would eventually strike the
         * plane.
         */
        bool InFrontOf(TVector3 direction);

    private:

        friend class GeometricTest;

        // A unit vector normal to the plane
        TVector3 normal;

        // The coefficient "d" in the equation ax + by + cz = d
        double coefficient;
    };

    /*
     * Represents a light ray with a time, position, and direction. All Ray objects travel at the speed of light in a
     * vacuum (units of cm/s).
     */
    class Ray
    {
    public:

        /*
         * The only constructor. Constructs a ray at the specified position, moving in the specified direction, at the
         * given time. If (0, 0, 0) is passed as the direction vector, (0, 0, 1) is used instead.
         */
        Ray(TVector3 position, TVector3 direction, double time);

        /*
         * Returns the current position of the ray.
         */
        TVector3 Position();

        /*
         * Returns the current velocity of the ray.
         */
        TVector3 Velocity();

        /*
         * Returns a unit vector pointing in the direction of the shower's motion.
         */
        TVector3 Direction();

        /*
         * When possible, use this method to set direction instead of directly modifying the velocity member.
         * This method ensures that the velocity vector has a magnitude equal to the speed of light in a vacuum.
         */
        void SetDirection(TVector3 direction);

        /*
         * Returns the current time of the ray.
         */
        double Time();

        /*
         * Moves the ray from its current position to the destination point. If the destination point doesn't lie along
         * the ray's current trajectory, the ray's direction is changed so it is pointing toward the destination.
         */
        void PropagateToPoint(TVector3 destination);

        /*
         * Moves the ray from its current position to the point where it intersects with the specified plane. If the ray
         * has already passed the plane, it will be moved BACKWARD. If the ray and the plane are exactly parallel, no
         * change is made to the ray's position.
         */
        void PropagateToPlane(Plane plane);

        /*
         * Find the point where this ray will intersect with the plane. If the ray and the plane are exactly parallel,
         * the current position of the ray is returned.
         */
        TVector3 PlaneImpact(Plane plane);

        /*
         * Finds the amount of time it will take before the ray will collide with the specified Plane object. Negative
         * times are returned if the ray has already passed the plane. Infinity is returned if the ray and plane are
         * exactly parallel.
         */
        double TimeToPlane(Plane plane);

        /*
         * Reflects the ray across the normal vector.
         */
        void Reflect(TVector3 normal);

        /*
         * Refracts the ray across the normal vector using the incident and outward indices of refraction specified. The
         * normal vector should point outward from the refracting surface, and therefore be opposite to the direction of
         * the ray's motion.
         */
        bool Refract(TVector3 normal, double n_in, double n_out);

        /*
         * Switches the the position and direction to a new, rotated reference frame.
         */
        void Transform(TRotation rotation);

    protected:

        friend class GeometricTest;

        // The current state of the Ray - cgs
        double time;
        TVector3 position;
        TVector3 velocity;

        /*
         * Moves the ray forward by the specified distance.
         */
        void IncrementPosition(double distance);

        /*
         * Move the time forward by some amount. This causes the ray to move forward according to its velocity.
         */
        void IncrementTime(double time_step);
    };

    /*
     * Represents an atmospheric cosmic ray shower. Contains various shower-specific parameters (x_0, n_max). Also
     * contains information about the density of the atmosphere. These atmospheric parameters are used to find the
     * relationship between slant depth and distance.
     */
    class Shower : public Ray
    {
    public:

        /*
         * A container for Shower parameters to be passed to the constructor.
         */
        struct Params
        {
            double energy;
            double x_max;
            double n_max;
            double rho_0;
            double scale_height;
            double delta_0;
        };

        /*
         * The only constructor. Takes all of the parameters used in the Ray class, as well as three shower-specific
         * parameters (x_0, x_max, n_max). Also takes two atmospheric parameters (rho_0 and scale_height). The starting
         * time is assumed to be zero by default.
         */
        Shower(Params params, TVector3 position, TVector3 direction, double time = 0);

        /*
         * Finds the age of the shower, given by 3 * X / (x + 2 * XMax).
         */
        double Age();

        /*
         * Returns the energy of the shower primary in MeV. Note that the energy is passed to the constructor in eV.
         */
        double EnergyMeV();

        /*
         * Returns the atmospheric density at the shower's current position.
         */
        double LocalRho();

        /*
         * Returns the local value for n - 1.
         */
        double LocalDelta();

        /*
         * Returns the number of particles in the shower at its current position.
         */
        double GaisserHillas();

        /*
         * Calculates the Cherenkov threshold energy of the shower.
         */
        double EThresh();

        /*
         * Increments the position of the shower by the specified slant depth. Returns the distance traversed.
         */
        void IncrementDepth(double depth);

    private:

        // Mask various methods from the Ray class. We have no need to change the shower's direction, and we only want
        // to change its position using the IncrementDepth method.
        using Ray::IncrementTime;
        using Ray::IncrementPosition;
        using Ray::SetDirection;
        using Ray::Reflect;
        using Ray::PropagateToPoint;
        using Ray::PropagateToPlane;

        friend class GeometricTest;

        // Parameters in the Gaisser-Hillas profile - cgs
        constexpr static double x_0 = -70.0;
        constexpr static double gh_lambda = 70.0;

        // Variable parameters in the Gaisser-Hillas profile - eV, cgs
        double energy;
        double x_max;
        double n_max;

        // Properties of the atmosphere - cgs
        double rho_0;
        double scale_height;
        double delta_0;

        /*
         * Returns the current depth of the shower in the atmosphere.
         */
        double X();
    };
}

#endif

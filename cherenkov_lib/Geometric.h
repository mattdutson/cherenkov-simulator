// Geometric.h
//
// Author: Matthew Dutson
//
// Defines the Plane, Shower, and Ray classes.

#ifndef GEOMETRIC_H
#define GEOMETRIC_H

#include <TRotation.h>
#include <TVector3.h>

#include "Utility.h"

namespace cherenkov_simulator
{
    /*
     * A class representing a plane in 3D space. Defined internally by a normal vector and the coefficient "d" in the 
     * equation ax + by + cz = d.
     */
    class Plane
    {
    public:

        /*
         * The default constructor. Sets the normal vector to (0, 0, 1) and the fixed point to the origin.
         */
        Plane();

        /*
         * Constructs a plane using the specified normal vector and intersection point. The point is converted to the 
         * coefficient "d" by plugging it into the equation ax + by + cz = d. Throws an invalid_argument exception if 
         * the normal vector is zero.
         */
        Plane(TVector3 normal, TVector3 point);

        /*
         * Returns a copy of the plane's normal vector.
         */
        TVector3 Normal() const;

        /*
         * Returns the coefficient "d" of the plane equation.
         */
        double Coefficient() const;

        /*
         * Returns true if a Ray going outward from the origin in the specified direction would eventually strike this
         * plane. If the plane is exactly at the origin, false is returned.
         */
        bool InFrontOf(TVector3 direction) const;

    private:

        friend class GeometricTest;

        TVector3 normal;
        double coeff;
    };

    /*
     * Represents a light Ray with a position, direction, and time. All Rays travel at the speed of light in a vacuum, 
     * measured in cm/s.
     */
    class Ray
    {
    public:

        /*
         * The default constructor. Sets the position to the origin, the direction to (0, 0, 1), and the time to zero.
         */
        Ray();

        /*
         * Constructs a Ray with the specified position, direction, and time. The direction does not necessarily have to
         * be unit (this will be taken care of by the constructor. An invalid_argument exception is thrown if the 
         * direction vector is zero.
         */
        Ray(TVector3 position, TVector3 direction, double cur_time);

        /*
         * Returns the current position of the Ray.
         */
        TVector3 Position() const;

        /*
         * Returns the current velocity of the Ray.
         */
        TVector3 Velocity() const;

        /*
         * Returns the unit vector of velocity.
         */
        TVector3 Direction() const;

        /*
         * Updates the direction with the one specified, setting the velocity to direction.Unit() * c. An 
         * invalid_argument exception is thrown if the direction vector is zero.
         */
        void SetDirection(TVector3 direction);

        /*
         * Returns the current time of the Ray.
         */
        double Time() const;

        /*
         * Moves the Ray from its current position to the destination point. If the destination doesn't lie along
         * the current trajectory, the Ray's direction is changed to the displacement between the current position and
         * the destination.
         */
        void PropagateToPoint(TVector3 destination);

        /*
         * Moves the Ray along its current trajectory until it collides with the Plane. If the Ray has already passed
         * the Plane, it will be moved backward, and its time will decrease. If the Ray and the Plane are exactly
         * parallel, no action is taken.
         */
        void PropagateToPlane(Plane plane);

        /*
         * Finds the point where this Ray will, or would have, collide with the Plane. If the Ray and the Plane are
         * exactly parallel, the current position of the Ray is returned.
         */
        TVector3 PlaneImpact(Plane plane) const;

        /*
         * Finds the amount of time it will take for the Ray to collide with the Plane. Negative times are returned if
         * the Ray has already passed the Plane. Infinity is returned if the Ray and Plane are exactly parallel.
         */
        double TimeToPlane(Plane plane) const;

        /*
         * Reflects the Ray across the the normal vector to some surface. The sign of the normal vector doesn't matter.
         * An invalid_argument exception is thrown if the normal vector is zero.
         */
        void Reflect(TVector3 normal);

        /*
         * Refracts the Ray across some interface with the normal vector, incident n, and outward n specified. The
         * normal vector should point outward from the refracting surface, opposite the direction of the ray. An
         * invalid_argument exception is thrown if the indices of refraction aren't at least 1 or if the normal vector
         * is zero. Returns true if the ray was refracted, and false if the ray was beyond the critical angle. If false
         * is returned, no changes are made to the vector's direction.
         */
        bool Refract(TVector3 normal, double n_in, double n_out);

        /*
         * Applies the rotation to both the Ray's position and velocity.
         */
        void Transform(TRotation rotation);

    protected:

        friend class GeometricTest;

        TVector3 position;
        TVector3 velocity;
        double cur_time;

        /*
         * Moves the Ray forward by the specified distance, updating the time in the process.
         */
        void IncrementPosition(double distance);

        /*
         * Move the time forward by some amount. This causes the Ray's position to change according to its velocity.
         */
        void IncrementTime(double time_step);
    };

    /*
     * Represents an atmospheric cosmic ray shower. Knows its own energy and the elevation of the detector (in order to
     * correctly calculate atmospheric densities).
     */
    class Shower : public Ray
    {
    public:

        /*
         * The default constructor. Objects constructed with this should only be used as placeholders.
         */
        Shower();

        /*
         * The main constructor. Takes the energy of the shower and the elevation of the detector. Also takes a
         * position, direction, and optional time, which are passed to the parent Ray constructor. An invalid_argument
         * exception is thrown if the energy is not positive.
         */
        Shower(double energy, double elevation, TVector3 position, TVector3 direction, double time = 0);

        /*
         * Finds the age of the shower, defined as 3 * X / (X + 2 * XMax).
         */
        double Age() const;

        /*
         * Returns the energy of the shower in MeV. Note that the energy is passed to the constructor in eV.
         */
        double EnergyMeV() const;

        /*
         * Returns the energy of the shower in eV.
         */
        double EnergyeV() const;

        /*
         * Calculates the impact parameter of the shower, assuming the detector is at the origin. See
         * http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html for an explanation of the point-line
         * distance formula.
         */
        double ImpactParam() const;

        /*
         * Calculates the angle of the shower with respect to a horizontal line in the shower-detector plane. Like
         * ImpactParam, this method assumes that the detector is at the origin.
         */
        double ImpactAngle() const;

        /*
         * Returns the atmospheric density at the shower's current position. The exponential model is defined at all
         * heights, so the function will return a value even if the height is negative.
         */
        double LocalRho() const;

        /*
         * Returns the local value of n - 1, where n is the index of refraction. This is proportional to the local
         * atmospheric density. For the same reason as LocalRho(), this function will return a value regardless of the
         * height of the shower.
         */
        double LocalDelta() const;

        /*
         * Uses the Gaisser-Hillas function to calculate the current number of particles (electrons) in the shower.
         */
        double GaisserHillas() const;

        /*
         * Calculates the Cherenkov threshold energy of the shower.
         */
        double EThresh() const;

        /*
         * Increases the slant depth of the shower by a specified amount, moving the shower forward in the process. The
         * assumption is made that the atmospheric density is constant over the course of the step.
         */
        void IncrementDepth(double depth);

        /*
         * Returns a header for rows of data created with ToString(). Used when writing CSV files.
         */
        static std::string Header();

        /*
         * Creates a string with a comma-separated impact parameter, impact angle, and shower direction. Used when
         * writing CSV files.
         */
        std::string ToString() const;

    private:

        // Mask unneeded Ray methods to keep them from being public.
        using Ray::IncrementTime;
        using Ray::IncrementPosition;
        using Ray::SetDirection;
        using Ray::Reflect;
        using Ray::PropagateToPoint;
        using Ray::PropagateToPlane;

        friend class GeometricTest;

        double energy;
        double elevation;

        /*
         * Calculates the shower's current slant depth.
         */
        double X() const;

        /*
         * Calculates the depth of the shower maximum.
         */
        double XMax() const;

        /*
         * Calculates the maximum number of particles in the shower.
         */
        double NMax() const;
    };
}

#endif

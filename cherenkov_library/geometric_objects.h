// geometric_objects.h
// cherenkov_library
//
// Created by Matthew Dutson on 9/8/16.
//
//

#ifndef geometric_objects_h
#define geometric_objects_h

#include "TVector3.h"

namespace cherenkov_simulator
{
    class Transformation
    {

    private:

        bool invert;

    public:

        /*
         * In the detector frame, the center of curvature is the origin. The axis of the telescope points in the +z
         * direction and the +x axis is parallel to the ground. The coordinate system is right-handed.
         */
        void TransformPosition(TVector3* position);

        void TransformDirection(TVector3* position);

        Transformation Inverse();
    };

    class Shower;

    class IntensityFunctor
    {
    public:

        virtual double GetIntensity(Shower shower) = 0;
    };

    class ConstantIntensity : public IntensityFunctor
    {
    public:
        double GetIntensity(Shower shower);
    };

    class Plane
    {
    private:

        TVector3 normal;

        double coefficient;

    public:

        Plane(TVector3 normal, TVector3 point);

        TVector3 Normal();

        double Coefficient();

        void Transform(Transformation t);
    };

    class Ray
    {
    protected:

        double current_time;

        TVector3 current_position;

        TVector3 current_velocity;

    public:

        Ray(double time, TVector3 position, TVector3 direction);

        void IncrementPosition(double time);

        double TimeToPlane(Plane p);

        double Time();

        TVector3 Position();

        TVector3 Velocity();

        void SetDirection(TVector3 direction);

        void Reflect(TVector3 normal);

        void PropagateToPoint(TVector3 point);

        void PropagateToPlane(Plane plane);

        void Transform(Transformation t);
    };

    class Shower : public Ray
    {
    private:

        double start_time;

        TVector3 start_position;

        IntensityFunctor intensity_functor;

    public:

        Shower(double time, TVector3 position, TVector3 direction, IntensityFunctor func);

        int NumberFluorescencePhotons();

        int NumberCherenkovPhotons();

        Ray GenerateCherenkovPhoton();

        void Transform(Transformation t);
    };
}

#endif

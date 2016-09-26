// geometric_objects.h
// cherenkov_simulator
//
// Created by Matthew Dutson on 9/8/16.
//
//

#ifndef geometric_objects_h
#define geometric_objects_h

#include "TVector3.h"

namespace cherenkov_simulator
{
    class Plane
    {
    private:
        
        TVector3 normal;
        
        double coefficient;
        
    public:
        
        Plane(TVector3 normal, TVector3 point);
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
    };

    class Shower: public Ray
    {
    private:
        
        double start_time;
        
        TVector3 start_position;
        
    public:
        
        Shower(double time, TVector3 position, TVector3 direction);
    };
}

#endif

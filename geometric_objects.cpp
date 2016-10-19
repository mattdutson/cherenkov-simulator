// geometric_objects.cpp
// cherenkov_simulator
//
// Created by Matthew Dutson on 9/8/16.
//
//

#include "geometric_objects.h"

namespace cherenkov_simulator
{

    double ConstantIntensity::GetIntensity(Shower shower) {
        return 0;
    }

    Ray::Ray(double time, TVector3 position, TVector3 direction) {

    }

    void Ray::IncrementPosition(double time) {

    }

    double Ray::TimeToPlane(Plane p) {
        return 0;
    }

    double Ray::GetTime() {
        return 0;
    }

    TVector3 Ray::GetPosition() {
        return TVector3();
    }

    TVector3 Ray::GetVelocity() {
        return TVector3();
    }

    void Ray::Reflect(TVector3 normal) {

    }

    void Ray::PropagateToPoint(TVector3 point) {

    }

    void Ray::PropagateToPlane(Plane plane) {

    }

    Shower::Shower(double time, TVector3 position, TVector3 direction, IntensityFunctor func): Ray(0, TVector3(), TVector3()) {
    }

    int Shower::NumberFluorescencePhotons() {
        return 0;
    }

    int Shower::NumberCherenkovPhotons() {
        return 0;
    }

    Ray Shower::GenerateCherenkovPhoton() {
        return Ray(0, TVector3(), TVector3());
    }

    Plane::Plane(TVector3 normal, TVector3 point) {

    }
}

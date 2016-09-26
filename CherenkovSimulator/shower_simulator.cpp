// shower_simulator.cpp
// cherenkov_simulator
//
// Created by Matthew Dutson on 9/8/16.
//
//

#include "simulator.h"
#include "geometric_objects.h"

namespace cherenkov_simulator
{
    SignalData Simulator::SimulateShower(Shower shower)
    {
        double time_delay = config.Get<double>("time_delay");
        Plane ground_plane = Plane(config.Get<TVector3>("ground_normal"), config.Get<TVector3>("ground_point"));
        int n_steps = (int) (shower.TimeToPlane(ground_plane) / time_delay) + 2;
        
        RawData raw_data = RawData();
        
        for(int i = 0; i < n_steps; i++) {
            ViewPoint(shower, &raw_data);
            shower.IncrementPosition(time_delay);
        }
        
        CherenkovBlast(shower, &raw_data);
        
        SignalData binned = BinRawData(raw_data);
        
        AddNoise(&binned);
        
        SignalData voltage = CalculateResponse(binned);
        
        return voltage;
    }
    
    void ViewPoint(Shower shower, RawData* raw_data)
    {
        TVector3 impact = RandomStopImpact();
        Ray observed = Ray(shower.GetTime(), shower.GetPosition(), impact - shower.GetPosition());
        
        Int_t numberOfPhotons = NumberOfPhotonsViewed(shower);
        for (Int_t i = 0; i < numberOfPhotons; i++) {
            TVector3 mirrorImpact = fMirror.GetMirrorImpact();
            TVector3 mirrorNormal = fMirror.GetMirrorNormal(mirrorImpact);
            TVector3 transformedPosition = shower.GetPosition();
            fCoordinates.PositionToObservatoryFrame(transformedPosition);
            TPlane3 focalPlane = TPlane3(TVector3(0, 0, 1), TVector3(0, 0, -fMirror.Radius() + FocalLength()));
            TRay detectedRay = TRay(shower.GetTime(), transformedPosition, mirrorImpact - transformedPosition);
            detectedRay.PropagateToPlane(focalPlane);
            if (CheckCollision(detectedRay.GetPosition())) {
                continue;
            }
            detectedRay.PropagateToPoint(mirrorImpact);
            detectedRay.ReflectFromPlane(TPlane3(mirrorNormal, mirrorImpact));
            detectedRay.PropagateToPlane(focalPlane);
            rawData.PushBack(detectedRay.GetPosition().X(), detectedRay.GetPosition().Y(), detectedRay.GetTime());
        }
    }
}

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
    
    void Simulator::ViewPoint(Shower shower, RawData* raw_data)
    {
        int total_photons = shower.GetPhotons();
        double fraction_captured = FractionCaptured(shower);
        int number_detected = (int) (total_photons * fraction_captured);
        
        for (int i = 0; i < number_detected; i++) {
            
            TVector3 impact = RandomStopImpact();
            Ray observed = Ray(shower.GetTime(), shower.GetPosition(), impact - shower.GetPosition());
            observed.PropagateToPoint(impact);
            
            observed.Refract(OuterLensNormal(impact), 1, config.Get<double>("refraction_index"));
            
            TVector3 exit;
            if (!LensExitPoint(observed, &exit)) continue;
            observed.PropagateToPoint(exit);
            
            observed.Refract(InnerLensNormal(exit), config.Get<double>("refraction_index"), 1);
            
            TVector3 reflect_point;
            if (!MirrorImpactPoint(observed, &reflect_point)) continue;
            if (CameraImpact(exit, reflect_point)) continue;
            observed.PropagateToPoint(reflect_point);
            
            observed.Reflect(MirrorNormal(reflect_point));
            
            TVector3 camera_impact;
            if (!CameraImpactPoint(observed, &camera_impact)) continue;
            observed.PropagateToPoint(camera_impact);
            
            TVector3 view_direction = GetViewDirection(camera_impact);
            
            raw_data->AddPoint(view_direction, observed.GetTime());
        }
    }
}

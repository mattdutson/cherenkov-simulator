// shower_simulator.cpp
// cherenkov_library
//
// Created by Matthew Dutson on 9/8/16.
//
//

#include "simulator.h"
#include "geometric_objects.h"
#include "TF1.h"
#include "TMath.h"

using namespace TMath;

namespace cherenkov_simulator
{
    Shower Simulator::GenerateRandomShower()
    {
        // Determine the energy of the shower primary using a distribution proportional to E^-energy_power.
        double min = config.Get("min_energy");
        double max = config.Get("max_energy");
        double power = config.Get("energy_power");
        std::string formula = "x^(-" + std::to_string(power) + ")";
        TF1 energy_profile = TF1("energy", formula.c_str(), min, max);
        double energy = energy_profile.GetRandom();

        // Determine the direction of the shower and its position relative to the detector. The angle of the shower
        // relative to the vertical goes as cos(theta) because shower have an isotropic flux in space.
        TF1 cosine_dist = TF1("cosine", "cos(x)", -TMath::Pi() / 2, TMath::Pi() / 2);
        TF1 uniform_angle_dist = TF1("uniform_angle", "1", 0, 2 * TMath::Pi());
        double theta = cosine_dist.GetRandom();
        double phi_shower = uniform_angle_dist.GetRandom();
        double phi_detector = uniform_angle_dist.GetRandom();

        // Determine the impact parameter.
        double r_min = config.Get("min_impact");
        double r_max = config.Get("max_impact");
        TF1 r_dist = TF1("r^2_dist", "x*x", r_min, r_max);
        double impact_param = r_dist.GetRandom();

        // Find the Cartesian shower axis vector. This vector is in the world frame (z is normal to the surface of the
        // earth, with x and y parallel to the surface. Note that the surface of the earth may not be parallel to the
        // local ground.
        TVector3 shower_axis = TVector3(cos(theta), sin(theta) * cos(phi_shower), sin(theta) * sin(phi_shower));

        // We define the origin of both the world and detector frames to be the detector's center of curvature for
        // simplicity.

        // We know that, at the impact point, the position vector of the shower is normal to its direction vector. Let's
        // find a vector normal to the direction. Choose the one with a zero z-component (we'll randomly rotate it
        // later). We solve dir_x * perp_x + dir_y * perp_y = 0.
        // TODO: Make sure this actually works
        TVector3 relative_direction;
        if (shower_axis.Y() == 0)
        {
            relative_direction = TVector3(0, 1, 0);
        } else
        {
            relative_direction = TVector3(1, shower_axis.X() / shower_axis.Y(), 0).Unit();
        }
        relative_direction.Rotate(phi_detector, shower_axis);

        // Now we have a point and direction for the shower (in the world frame). We will now determine other shower
        // parameters, assuming a proton primary.
        // TODO: Add specific numerical values to the configuration file.

        // Find the depth of the first interaction. See AbuZayyad 6.1
        TF1 depth_dist = TF1("interact", "e^(-x/70)", 0, TMath::Infinity());
        double first_interaction = depth_dist.GetRandom();

        // Determine the depth of the shower maximum. See AbuZayyad 6.2
        double x_max_diff = 725.0 + 55.0 * (Log(energy) - 18.0) - 70;

        // Determine the size of the shower at the maximum. See AbuZayyad 6.4
        double n_max = energy / (1.3e9);
    }

    VoltageSignal Simulator::SimulateShower(Shower shower)
    {
        // Transform the shower to the detector frame.
        shower.Transform(to_detector_frame);

        // Construct the ground plane and transform it to the detector reference frame.
        Plane ground_plane = Plane(config.Get<TVector3>("ground_normal"), config.Get<TVector3>("ground_point"));
        ground_plane.Transform(to_detector_frame);

        // Determine the number of steps in the simulation of the shower
        double time_delay = config.Get<double>("time_delay");
        int n_steps = (int) (shower.TimeToPlane(ground_plane) / time_delay) + 2;

        // Initialize the data container for photon counts.
        PhotonCount photon_count = PhotonCount(config);
        
        for(int i = 0; i < n_steps; i++) {
            ViewFluorescencePhotons(shower, &photon_count);
            ViewCherenkovPhotons(shower, ground_plane, &photon_count);
            shower.IncrementPosition(time_delay);
        }
        
        AddNoise(&photon_count);
        
        return VoltageResponse(photon_count);
    }
    
    void Simulator::ViewFluorescencePhotons(Shower shower, PhotonCount* photon_count)
    {
        int total_photons = shower.NumberFluorescencePhotons();
        double fraction_captured = FluorescenceFractionCaptured(shower);
        int number_detected = (int) (total_photons * fraction_captured);
        
        for (int i = 0; i < number_detected; i++)
        {
            TVector3 lens_impact = RandomStopImpact();
            Ray photon = Ray(shower.Time(), shower.Position(), lens_impact - shower.Position());
            photon.PropagateToPoint(lens_impact);
            
            SimulateOptics(photon, photon_count);
        }
    }
    
    void Simulator::ViewCherenkovPhotons(Shower shower, Plane ground_plane, PhotonCount* photon_count)
    {
        int cherenkov_photons = shower.NumberCherenkovPhotons();
        
        for (int i = 0; i < cherenkov_photons; i++)
        {
            Ray photon = shower.GenerateCherenkovPhoton();
            photon.PropagateToPlane(ground_plane);
            if (!ReflectFromGround(&photon)) continue;
            TVector3 stop_impact;
            if (!LensImpactPoint(photon, &stop_impact)) continue;
            photon.PropagateToPoint(stop_impact);
            
            SimulateOptics(photon, photon_count);
        }
    }
    
    void Simulator::SimulateOptics(Ray photon, PhotonCount* photon_count)
    {
        DeflectFromLens(&photon);
        
        TVector3 reflect_point;
        if (!MirrorImpactPoint(photon, &reflect_point)) return;
        if (BlockedByCamera(photon.Position(), reflect_point)) return;
        photon.PropagateToPoint(reflect_point);
        photon.Reflect(MirrorNormal(reflect_point));
        
        TVector3 camera_impact;
        if (!CameraImpactPoint(photon, &camera_impact)) return;
        photon.PropagateToPoint(camera_impact);
        
        photon_count->AddPoint(GetViewDirection(camera_impact), photon.Time());
    }

    void Simulator::AddNoise(PhotonCount *photon_count) {

    }

    VoltageSignal Simulator::VoltageResponse(PhotonCount photon_count) {
        return VoltageSignal();
    }

    TVector3 Simulator::RandomStopImpact() {
        return TVector3();
    }

    double Simulator::FluorescenceFractionCaptured(Shower shower) {
        return 0;
    }

    TVector3 Simulator::MirrorNormal(TVector3 point) {
        return TVector3();
    }

    bool Simulator::LensImpactPoint(Ray ray, TVector3 *point) {
        return false;
    }

    bool Simulator::MirrorImpactPoint(Ray ray, TVector3 *point) {
        return false;
    }

    bool Simulator::CameraImpactPoint(Ray ray, TVector3 *point) {
        return false;
    }

    bool Simulator::BlockedByCamera(TVector3 start, TVector3 end) {
        return false;
    }

    void Simulator::DeflectFromLens(Ray *photon) {

    }

    bool Simulator::ReflectFromGround(Ray *photon) {
        return false;
    }

    void Simulator::ImpactPointToCameraIndex(TVector3 impact, int *x_index, int *y_index) {

    }

    TVector3 Simulator::GetViewDirection(TVector3 impact_point) {
        return TVector3();
    }

    Shower Simulator::ReconstructShower(VoltageSignal dat) {
        return Shower(0, TVector3(), TVector3(), ConstantIntensity());
    }

    void Simulator::EstimateAccuracy() {

    }
}

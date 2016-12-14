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
#include "TRandom3.h"
#include "common.h"

using namespace TMath;

namespace cherenkov_simulator
{
    Simulator::Simulator(FileOptions* file_config)
    {
        // The overall program configuration
        config = file_config;

        // A general-purpose random number generator
        rng = TRandom3();

        // The rotation from detector frame to world frame. The frames share the same x-axis.
        rotate_to_world = TRotation();
        rotate_to_world.RotateX(-PiOver2() + config->Get<double>("elevation_angle"));

        // The distribution of shower energies
        std::string energy_formula = "x^(-" + std::to_string(config->Get<double>("energy_power")) + ")";
        energy_distribution = TF1("energy", energy_formula.c_str(), config->Get<double>("min_energy"),
                                  config->Get<double>("max_energy"));

        // The distribution of shower vertical directions
        cosine_distribution = TF1("cosine", "cos(x)", -TMath::Pi() / 2, TMath::Pi() / 2);

        // The Monte Carlo distribution of impact parameters
        impact_distrubition = TF1("impact", "x^2", config->Get<double>("min_impact"),
                                  config->Get<double>("max_impact"));

        // The distribution of first interaction depths (See AbuZayyad 6.1)
        std::string interaction_formula = "e^(-x/ " + std::to_string(config->Get<double>("mean_interaction")) + ")";
        interact_distribution = TF1("interact", interaction_formula.c_str(), 0, TMath::Infinity());

        // Atmospheric parameters (see notes)
        H = config->Get<double>("atmosphere.H");
        double sea_density = config->Get<double>("atmosphere.sea_level");
        double elev = config->Get<double>("atmosphere.detector_elevation");
        rho_0 = sea_density * Exp(-elev / H);
    }

    VoltageSignal Simulator::SimulateShower(Shower shower)
    {
        // Construct the ground plane in the world frame.
        Plane ground_plane = Plane(config->Get<TVector3>("ground_normal"), config->Get<TVector3>("ground_point"));

        // Initialize the data container for photon counts.
        PhotonCount photon_count = PhotonCount(config);

        // Step the shower through its path.
        int n_steps = config->Get<int>("simulation.n_steps");
        for(int i = 0; i < n_steps; i++) {
            ViewFluorescencePhotons(shower, &photon_count);
            ViewCherenkovPhotons(shower, ground_plane, &photon_count);
            double distance = FindDistance(shower);
            shower.IncrementPosition(distance);
        }
        
        AddNoise(&photon_count);
        
        return VoltageResponse(photon_count);
    }

    double Simulator::FindDistance(Shower shower)
    {
        double total_depth = VerticalDepth(shower.GroundImpact(), shower.StartPosition());
        double depth_step = total_depth / n_steps;
        double vertical_distance = -H * Log(Exp(-shower.Position().Z() / H) + H * depth_step / rho_0);
        return vertical_distance / Abs(shower.Velocity().CosTheta());
    }

    Shower Simulator::GenerateRandomShower()
    {
        double energy = energy_distribution.GetRandom();

        // Determine the direction of the shower and its position relative to the detector. The angle of the shower
        // relative to the vertical goes as cos(theta) because shower have an isotropic flux in space.
        double theta = cosine_distribution.GetRandom();
        double phi_shower = rng.Uniform(2 * Pi());

        // Determine the impact parameter.
        double impact_param = impact_distrubition.GetRandom();

        // Find the Cartesian shower axis vector. This vector is in the world frame (z is normal to the surface of the
        // earth, with x and y parallel to the surface. Note that the surface of the earth may not be parallel to the
        // local ground.
        TVector3 shower_axis = TVector3(cos(theta), sin(theta) * cos(phi_shower), sin(theta) * sin(phi_shower));

        // We define the origin of both the world and detector frames to be the detector's center of curvature for
        // simplicity. We know that, at the impact point, the position vector of the shower is normal to its direction
        // vector. Let's find a random vector normal to the shower axis.
        TVector3 impact_point = impact_param * RandomPerpendicularVector(shower_axis, rng);

        // Find the depth of the first interaction, the depth of the maximum, and the size of the shower maximum (See
        // AbuZayyad 6.1-6.4). We assume a proton primary.
        // TODO: Add specific numerical values to the configuration file.
        double x_0 = interact_distribution.GetRandom();
        double x_max = 725.0 + 55.0 * (Log(energy) - 18.0) - 70 + x_0;
        double n_max = energy / config->Get<double>("n_max_ratio");

        // Trace the path of the shower back to the location of the first interaction. Start by finding the elevation of
        // the first interaction.
        double interaction_height = -H * cos(theta) * Log(x_0 / (rho_0 * H * cos(theta)));
        double param = (interaction_height - impact_point.Z()) / (shower_axis.Z());
        TVector3 starting_position = impact_point + param * shower_axis;

        // Create a new shower with all of the randomly determined parameters.
        return Shower(starting_position, shower_axis, x_0, x_max, n_max);
    }
    
    void Simulator::ViewFluorescencePhotons(Shower shower, PhotonCount* photon_count)
    {
        int number_detected = NumberFluorescencePhotons(shower);
        
        for (int i = 0; i < number_detected; i++)
        {
            TVector3 lens_impact = RandomStopImpact();
            Ray photon = Ray(shower.Time(), shower.Position(), lens_impact - shower.Position());
            photon.PropagateToPoint(lens_impact);
            
            SimulateOptics(photon, photon_count);
        }
    }

    int Simulator::NumberFluorescencePhotons(Shower shower)
    {
        return (int) GaiserHilles(shower) * config->Get<double>("fluorescence_yield") *
               PhotonFraction(shower.Position());
    }

    int Simulator::NumberCherenkovPhotons(Shower shower)
    {
        // TODO: Use the Gaiser-Hilles profile and material from Nerling to determine the number of Cherenkov photons
        int total_produced = 0;

        // Determine the fraction captured.
        return total_produced * PhotonFraction(shower.GroundImpact());
    }

    double Simulator::GaiserHilles(Shower shower)
    {
        // Use the Gaisser-Hilles profile to find the number of electrons.
        double x = SlantDepth(shower.Position(), shower.StartPosition());
        double n_max = shower.NMax();
        double x_max = shower.XMax();
        double x_0 = shower.X0();
        return n_max * Power((x - x_0) / (x_max - x_0), (x_max - x_0) / lambda) * Exp((x_max - x) / lambda);
    }

    double Simulator::PhotonFraction(TVector3 view_point)
    {
        TVector3 detector_axis = rotate_to_world * TVector3(0, 0, 1);
        double angle = detector_axis.Angle(view_point);
        double distance = view_point.Mag();
        double area_fraction = stop_radius * stop_radius / (4.0 * distance * distance);
        return area_fraction * cos(angle);
    }

    Ray Simulator::GenerateCherenkovPhoton(Shower shower)
    {
        TVector3 rotation_axis = RandomPerpendicularVector(shower.Velocity().Unit(), rng);
        std::string formula = "e^(-x/" + std::to_string(theta_0) + ")/sin(x)";
        TF1 angular_distribution = TF1("distro", formula.c_str(), 0, Pi());
    }

    double Simulator::VerticalDepth(TVector3 point1, TVector3 point2)
    {
        if (point1.Z() > point2.Z())
        {
            // point1 should be the point with the smaller z-coordinate
            TVector3 tmp = point1;
            point1 = point2;
            point2 = tmp;
        }
        return -rho_0 / H * (Exp(-point2.Z() / H) - Exp(-point1.Z() / H));
    }

    double Simulator::SlantDepth(TVector3 point1, TVector3 point2)
    {
        TVector3 diff = point1 - point2;
        double cosine = Abs(diff.Z() / diff.Mag());
        return VerticalDepth(point1, point2) / cosine;
    }
    
    void Simulator::ViewCherenkovPhotons(Shower shower, Plane ground_plane, PhotonCount* photon_count)
    {
        int number_detected = NumberCherenkovPhotons(shower);

        for (int i = 0; i < number_detected; i++)
        {
            Ray photon = GenerateCherenkovPhoton(shower);
            photon.PropagateToPlane(ground_plane);
            TVector3 stop_impact;
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
        return Shower(TVector3(), TVector3(), 0, 0, 0);
    }

    void Simulator::EstimateAccuracy() {

    }
}

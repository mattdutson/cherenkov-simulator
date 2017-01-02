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
#include "Math/Polynomial.h"

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
        double phi_shower = rng.Uniform(TwoPi());

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
        // Refract the ray over the corrector plate.
        DeflectFromLens(&photon);

        // Check whether the ray hit the back of the camera.
        TVector3 reflect_point, camera_impact;
        if (CameraImpactPoint(photon, &camera_impact)) return;

        // Check whether the ray bounces off the mirror and then reflect it.
        if (!MirrorImpactPoint(photon, &reflect_point)) return;
        photon.PropagateToPoint(reflect_point);
        photon.Reflect(MirrorNormal(reflect_point));

        // Check whether and where the ray is detected by the cluster.
        if (!CameraImpactPoint(photon, &camera_impact)) return;
        photon.PropagateToPoint(camera_impact);

        // Record the detection position.
        photon_count->AddPoint(GetViewDirection(camera_impact), photon.Time());
    }

    void Simulator::AddNoise(PhotonCount *photon_count) {

    }

    VoltageSignal Simulator::VoltageResponse(PhotonCount photon_count) {
        return VoltageSignal();
    }

    TVector3 Simulator::RandomStopImpact() {
        double phi_rand = rng.Uniform(TwoPi());

        // The probability that a random point on a circle will lie at radius r is proportional to r (area of slice is
        // 2*pi*r*dr).
        TF1 radius_distribution = TF1("stop_impact", "x", 0, config->Get<double>("camera.stop_size"));
        double r_rand = radius_distribution.GetRandom();
        return TVector3(r_rand * Cos(phi_rand), r_rand * Sin(phi_rand), 0);
    }

    TVector3 Simulator::MirrorNormal(TVector3 point) {
        // For a spherical mirror, the normal vector will always point straight back to the center of curvature (the
        // origin in this case).
        return -point;
    }

    bool Simulator::LensImpactPoint(Ray ray, TVector3 *point) {
        return false;
    }

    bool Simulator::MirrorImpactPoint(Ray ray, TVector3* point)
    {
        // Find the point with the smallest z-component where the ray intersects with the mirror sphere
        double radius = config->Get<double>("camera.mirror_radius");
        BackOriginSphereImpact(ray, point, radius);

        // Ensure that the ray actually hits the mirror
        double mirror_size = config->Get("camera.mirror_size") / 2.0;
        return WithinXYDisk(*point, mirror_size) && point->Z() < 0;
    }

    bool Simulator::BackOriginSphereImpact(Ray ray, TVector3* point, double radius)
    {
        // Solve for the time when the ray impacts a sphere centered at the origin. See notes for details.
        double a = ray.Position().Dot(ray.Position()) - radius * radius;
        double b = 2 * ray.Position().Dot(ray.Velocity());
        double c = ray.Velocity().Dot(ray.Velocity());
        ROOT::Math::Polynomial poly = ROOT::Math::Polynomial(a, b, c);
        std::vector<double> roots = poly.FindRealRoots();
        if (roots.size() == 0)
        {
            *point = TVector3(0, 0, 0);
            return false;
        }
        double time;
        if (roots.size() == 1)
        {
            time = roots[0];
        } else
        {
            if (ray.Velocity().Z() < 0)
            {
                time = Max(roots[0], roots[1]);
            } else
            {
                time = Min(roots[0], roots[1]);
            }
        }
        *point = ray.Position() + time * ray.Velocity();
        return true;
    }

    bool Simulator::WithinXYDisk(TVector3 vec, double radius)
    {
        double xy_radius = Sqrt(vec.X() * vec.X() + vec.Y() * vec.Y());
        if (xy_radius < radius)
        {
            return true;
        } else
        {
            return false;
        }
    }

    bool Simulator::CameraImpactPoint(Ray ray, TVector3 *point) {
        double radius = config->Get<double>("camera.mirror_radius") / 2.0;
        BackOriginSphereImpact(ray, point, radius);

        double cluster_size = config->Get<double>("camera.cluster_size") / 2.0;
        return WithinXYDisk(*point, cluster_size) && point->Z() < 0;
    }

    void Simulator::DeflectFromLens(Ray *photon) {
        double stop_radius = config->Get<double>("stop_size");
        TVector3 position = photon->Position();

        // Only deflect the photon if it's on the portion of the lens where the s^4 term dominates.
        double photon_axis_dist = Sqrt(position.X() * position.X() + position.Y() * position.Y());
        if (photon_axis_dist < Sqrt2() * stop_radius)
        {
            return;
        }

        double refrac_corrector = config->Get<double>("optics.n_corrector");
        double refrac_air = config->Get<double>("optics.n_air");
        double r_curve = config->Get<double>("camera.mirror_radius");
        double schmidt_coefficient = 1.0 / (4.0 * (refrac_corrector - 1) * Power(r_curve, 3));

        // Now let's find a vector normal to the lens surface. Start with the derivative.
        double deriv = 4 * schmidt_coefficient * Power(photon_axis_dist, 3);
        double theta = ATan(1.0 / deriv);

        double phi = position.Phi();

        TVector3 norm = TVector3(Sin(theta) * Cos(phi), Sin(theta) * Sin(phi), Cos(phi));

        photon->Refract(norm, refrac_air, refrac_corrector);

        photon->Refract(TVector3(0, 0, 1), refrac_corrector, refrac_air);

    }

    bool Simulator::ReflectFromGround(Ray *photon) {
        return false;
    }

    void Simulator::ImpactPointToCameraIndex(TVector3 impact, int *x_index, int *y_index) {
        TVector3 direction = GetViewDirection(impact);

        // The angle of the yz projection with the z axis. See notes for details.
        double elevation = ATan(direction.Y() / direction.Z());

        // The angle of the xz projection with the z axis. See notes for details.
        double azimuth = ATan(direction.X() / direction.Z());

        // Arc length = angle * radius
        double cluster_size = config->Get<double>("cluster_size");
        double mirror_radius = config->Get<double>("mirror_radius");
        double field_of_view = cluster_size / (mirror_radius / 2.0);

        int n = config->Get<int>("n_pmt_across");

        double y_bin_size = field_of_view / ((double) n);
        *y_index = Floor(elevation / y_bin_size) + n + 1;

        double x_bin_size = field_of_view / ((double) n * Cos(elevation));
        *x_index = Floor(azimuth / x_bin_size) + n + 1;
    }

    TVector3 Simulator::GetViewDirection(TVector3 impact_point) {
        // TODO: Investigate the effect of the corrector plate on the outward view direction
        return -impact_point.Unit();
    }

    Shower Simulator::ReconstructShower(VoltageSignal dat) {
        return Shower(TVector3(), TVector3(), 0, 0, 0);
    }

    void Simulator::EstimateAccuracy() {

    }
}

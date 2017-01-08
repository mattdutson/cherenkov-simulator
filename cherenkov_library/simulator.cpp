// shower_simulator.cpp
// cherenkov_library
//
// Created by Matthew Dutson on 9/8/16.
//
//

#include "simulator.h"
#include "Math/Polynomial.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

using namespace TMath;
using std::string;

namespace cherenkov_simulator
{
    Simulator::Simulator() : ground_plane(TVector3(0, 0, 0), TVector3(0, 0, 0))
    {
    }

    void Simulator::ParseFile(std::ifstream config_file)
    {
        using namespace boost::property_tree;

        // Parse the file to XML.
        ptree xml_file = ptree();
        try
        {
            read_xml(config_file, xml_file);
        }
        catch (...)
        {
            throw std::runtime_error("There was a problem parsing the file to XML. Check for syntax errors.");
        }
        ptree config = xml_file.get_child("configuration");

        // Parameters related to the behavior of the simulation
        {
            depth_step = config.get<double>("depth_step");
        }

        // Parameters relating to the position and orientation of the detector relative to its surroundings
        {
            // Construct the ground plane in the world frame.
            ground_plane = Plane(ToVector(config.get<string>("ground_normal")),
                                 ToVector(config.get<string>("ground_point")));

            // The rotation from detector frame to world frame. The frames share the same x-axis.
            rotate_to_world = TRotation();
            rotate_to_world.RotateX(-PiOver2() + config.get<double>("elevation_angle"));
        }

        // Parameters used to generate random showers in the Monte Carlo simulation
        {
            // The distribution of shower energies
            std::string energy_formula = "x^(-" + std::to_string(config.get<double>("energy_pow")) + ")";
            energy_distribution = TF1("energy", energy_formula.c_str(), config.get<double>("e_min"),
                                      config.get<double>("e_max"));

            // The distribution of shower vertical directions
            cosine_distribution = TF1("cosine", "cos(x)", -TMath::Pi() / 2, TMath::Pi() / 2);

            // The Monte Carlo distribution of impact parameters
            impact_distribution = TF1("impact", "x", config.get<double>("impact_min"),
                                      config.get<double>("impact_max"));

            // The distribution of first interaction depths (See AbuZayyad 6.1)
            std::string interaction_formula = "e^(-x/ " + std::to_string(config.get<double>("avg_interact")) + ")";
            interact_distribution = TF1("interact", interaction_formula.c_str(), 0, TMath::Infinity());
        }

        // Parameters defining properties of the atmosphere
        {
            scale_height = config.get<double>("atmosphere.scale_height");
            double rho_sea = config.get<double>("atmosphere.rho_sea");
            double detect_elevation = config.get<double>("atmosphere.detect_elevation");
            rho_0 = rho_sea * Exp(-detect_elevation / scale_height);

            // 1 - the index of refraction at the detector (proportional to atmospheric density)
            double delta_sea = 1.0 - config.get<double>("optics.n_air");
            delta_0 = delta_sea * Exp(-detect_elevation / scale_height);
        }

        // Physics constants
        {
            mass_e = config.get<double>("mass_e");
            fine_struct = config.get<double>("fine_struct");
        }

        // A general-purpose random number generator
        rng = TRandom3();
    }

    double Simulator::AtmosphereDensity(double height)
    {
        return rho_0 * Exp(-height / scale_height);
    }

    double Simulator::AtmosphereDelta(double height)
    {
        // Assume that delta is proportional to the local atmospheric density.
        return delta_0 * Exp(-height / scale_height);
    }

    double Simulator::EThresh(double height)
    {
        double delta = AtmosphereDelta(height);
        return mass_e * C() * C() / Sqrt(2 * delta);
    }

    double Simulator::ThetaC(double height)
    {
        return ckv_k1 * Power(EThresh(height), ckv_k2);
    }

    VoltageSignal Simulator::SimulateShower(Shower shower)
    {
        // TODO: Change the initial time in the data container to a better value.
        // Arc length = angle * radius
        double pmt_angular_size = (cluster_size / (double) n_pmt_across) / (mirror_radius / 2.0);
        PhotonCount photon_count = PhotonCount(n_pmt_across, shower.Time(), time_bin, pmt_angular_size);

        // Step the shower through its path.
        while (shower.TimeToPlane(ground_plane) > 0)
        {
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
        double vertical_distance =
                -scale_height * Log(Exp(-shower.Position().Z() / scale_height) + scale_height * depth_step / rho_0);
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
        double impact_param = impact_distribution.GetRandom();

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
        double n_max = energy / n_max_ratio;

        // Trace the path of the shower back to the location of the first interaction. Start by finding the elevation of
        // the first interaction.
        double interaction_height = -scale_height * cos(theta) * Log(x_0 / (rho_0 * scale_height * cos(theta)));
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
        int n_charged = GaiserHilles(shower);
        double alpha_eff = IonizationLossRate(shower);
        double yield = FluorescenceYield(shower);

        // See Stratton Equation 4.2. The energy deposit rate for a single photon is alpha_eff, so the deposit rate for
        // all photons is alpha_eff * N.
        double total_produced = n_charged * alpha_eff * yield * depth_step;

        // Find the fraction captured by the camera.
        return total_produced * PhotonFraction(shower.Position());
    }

    int Simulator::NumberCherenkovPhotons(Shower shower)
    {
        // 2 * pi * alpha * (1 / lambda1 - 1 / lambda2) / rho
        double rho = AtmosphereDensity(shower.Position().Z());
        double k_out = 2 * Pi() * fine_struct / rho * (1 / lambda_min - 1 / lambda_max);
        double k_1 = k_out * 2 * AtmosphereDelta(shower.Position().Y());

        // TODO: Make sure the value for C is in units of cm/s.
        double k_2 = k_out * Sq(mass_e) * Power(C(), 4);

        // Parameters in the electron energy distribution
        double age = ShowerAge(shower);
        double a1 = fe_a11 - fe_a12 * age;
        double a2 = fe_a21 - fe_a22 * age;
        double a0 = fe_k0 * Exp(fe_k1 * age + fe_k2 * Sq(age));

        // See notes for details. Integrate over lnE from lnE_thresh to infinity. E terms turn to e^E because the
        // variable of integration is lnE.
        std::stringstream func_string;
        func_string << "(" << a0 << "* e^x / ((" << a1 << "+e^x)(" << a2 << "+e^x)^" << age << "))(" << k_1 << "-"
                    << k_2 << "/e^x)";
        double e_thresh = EThresh(shower.Position().Z());
        TF1 func = TF1("integrand", func_string.str().c_str(), e_thresh, Infinity());
        double integral = func.Integral(e_thresh, Infinity());

        // Determine the number captured. Multiply the PhotonFraction by two because we're dealing with a half sphere,
        // not a full sphere.
        return GaiserHilles(shower) * integral * PhotonFraction(shower.GroundImpact()) * 2;
    }

    double Simulator::GaiserHilles(Shower shower)
    {
        // Use the Gaisser-Hilles profile to find the number of electrons.
        double x = SlantDepth(shower.Position(), shower.StartPosition());
        double n_max = shower.NMax();
        double x_max = shower.XMax();
        double x_0 = shower.X0();
        return n_max * Power((x - x_0) / (x_max - x_0), (x_max - x_0) / gh_lambda) * Exp((x_max - x) / gh_lambda);
    }

    double Simulator::FluorescenceYield(Shower shower)
    {
        // TODO: Implement fluorescence yield calculation
        return 0.0;
    }

    double Simulator::ShowerAge(Shower shower)
    {
        double x = SlantDepth(shower.Position(), shower.StartPosition());
        return 3.0 * x / (x + 2.0 * shower.XMax());
    }

    double Simulator::IonizationLossRate(Shower shower)
    {
        double s = ShowerAge(shower);
        return ion_c1 / Power(ion_c2 + s, ion_c3) + ion_c4 + ion_c5 * s;
    }

    double Simulator::PhotonFraction(TVector3 view_point)
    {
        TVector3 detector_axis = rotate_to_world * TVector3(0, 0, 1);
        double angle = detector_axis.Angle(view_point);
        double distance = view_point.Mag();
        double area_fraction = Sq(stop_size / 2.0) / (4.0 * Sq(distance));
        return area_fraction * cos(angle);
    }

    Ray Simulator::GenerateCherenkovPhoton(Shower shower)
    {
        TVector3 rotation_axis = RandomPerpendicularVector(shower.Velocity().Unit(), rng);
        std::string formula = "e^(-x/" + std::to_string(ThetaC(shower.Position().Z())) + ")/sin(x)";
        TF1 angular_distribution = TF1("distro", formula.c_str(), 0, Pi());
        TVector3 direction = shower.Velocity().Unit();
        direction.Rotate(angular_distribution.GetRandom(), rotation_axis);
        return Ray(shower.Time(), shower.Position(), direction);
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
        return -rho_0 / scale_height * (Exp(-point2.Z() / scale_height) - Exp(-point1.Z() / scale_height));
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
        photon_count->AddPhoton(photon.Time(), camera_impact);
    }

    void Simulator::AddNoise(PhotonCount* photon_count)
    {
        // TODO: Find expressions for these in terms of the detector field of view and the radius of curvature
        double detector_area;
        double detector_solid_angle;
        double pmt_area;
        double pmt_solid_angle;
    }

    VoltageSignal Simulator::VoltageResponse(PhotonCount photon_count) {
        return VoltageSignal();
    }

    TVector3 Simulator::RandomStopImpact() {
        double phi_rand = rng.Uniform(TwoPi());

        // The probability that a random point on a circle will lie at radius r is proportional to r (area of slice is
        // 2*pi*r*dr).
        TF1 radius_distribution = TF1("stop_impact", "x", 0, stop_size / 2.0);
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
        BackOriginSphereImpact(ray, point, mirror_radius);

        // Ensure that the ray actually hits the mirror
        return WithinXYDisk(*point, mirror_size / 2.0) && point->Z() < 0.0;
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
        BackOriginSphereImpact(ray, point, mirror_radius);
        return WithinXYDisk(*point, cluster_size / 2.0) && point->Z() < 0;
    }

    void Simulator::DeflectFromLens(Ray *photon) {
        TVector3 position = photon->Position();

        // Only deflect the photon if it's on the portion of the lens where the s^4 term dominates.
        double photon_axis_dist = Sqrt(position.X() * position.X() + position.Y() * position.Y());
        if (photon_axis_dist < Sqrt2() * stop_size / 2.0)
        {
            return;
        }

        double schmidt_coefficient = 1.0 / (4.0 * (refrac_lens - 1) * Power(mirror_radius, 3));

        // Now let's find a vector normal to the lens surface. Start with the derivative.
        double deriv = 4 * schmidt_coefficient * Power(photon_axis_dist, 3);
        double theta = ATan(1.0 / deriv);

        double phi = position.Phi();

        TVector3 norm = TVector3(Sin(theta) * Cos(phi), Sin(theta) * Sin(phi), Cos(phi));

        photon->Refract(norm, delta_0 + 1, refrac_lens);

        photon->Refract(TVector3(0, 0, 1), refrac_lens, delta_0 + 1);

    }

    bool Simulator::ReflectFromGround(Ray *photon) {
        return false;
    }

    // TODO: Move this method to the photon count data container.
    void Simulator::ImpactPointToCameraIndex(TVector3 impact, int *x_index, int *y_index) {

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

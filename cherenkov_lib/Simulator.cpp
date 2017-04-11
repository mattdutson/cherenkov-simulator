// Simulator.cpp
//
// Author: Matthew Dutson
//
// Implementation of Simulator.h

#include <boost/property_tree/ptree.hpp>
#include <Math/Polynomial.h>
#include <TMath.h>

#include "Simulator.h"
#include "Utility.h"

using namespace TMath;

using std::string;
using boost::property_tree::ptree;

namespace cherenkov_simulator
{
    Simulator::Simulator(const ptree& config)
    {
        // Behavior of the simulation
        depth_step = config.get<double>("simulation.depth_step");
        fluor_thin = config.get<double>("simulation.fluor_thin");
        ckv_thin = config.get<double>("simulation.ckv_thin");
        back_toler = config.get<double>("simulation.back_toler");

        // Surrounding setup and orientation
        TVector3 ground_normal = Utility::ToVector(config.get<string>("surroundings.ground_normal"));
        TVector3 ground_point = Utility::ToVector(config.get<string>("surroundings.ground_point"));
        ground_plane = Plane(ground_normal, ground_point);
        rotate_to_world = Utility::MakeRotation(config.get<double>("surroundings.elevation_angle"));

        // Properties of the detector
        mirror_radius = config.get<double>("detector.mirror_radius");
        stop_diameter = mirror_radius / (2.0 * config.get<double>("detector.f_number"));
        double field_of_view = config.get<double>("detector.field_of_view");
        mirror_size = stop_diameter + 2.0 * mirror_radius * Tan(field_of_view / 2.0);
        cluster_diameter = mirror_radius * Sin(field_of_view / 2.0);

        // PhotonCount parameters
        count_params.bin_size = config.get<double>("simulation.time_bin");
        count_params.n_pixels = config.get<int>("detector.n_pmt_across");
        count_params.linear_size = cluster_diameter / count_params.n_pixels;
        count_params.angular_size = count_params.linear_size / (mirror_radius / 2.0);

        // Miscellaneous member variables
        ckv_integrator = TF1("ckv_integrator", ckv_func, 0.0, Infinity(), 3);
        ckv_integrator.SetParNames("age", "rho", "delta");
        rng = TRandom3();
        if (config.get<bool>("simulation.time_seed")) rng.SetSeed();
    }

    PhotonCount Simulator::SimulateShower(Shower shower)
    {
        // A lower bound on the time the first photon will reach the detector
        double time = shower.Time() + (shower.Position().Mag() - mirror_radius) / Utility::c_cent;
        PhotonCount photon_count = PhotonCount(count_params, MinTime(shower), MaxTime(shower));

        // Step the shower along its path
        while (shower.TimeToPlane(ground_plane) > 0)
        {
            shower.IncrementDepth(depth_step);
            ViewFluorescencePhotons(shower, photon_count);
            ViewCherenkovPhotons(shower, ground_plane, photon_count);
        }
        photon_count.Trim();
        return photon_count;
    }

    double Simulator::CherenkovFunc::operator()(double* x, double* p)
    {
        // See Simulator constructor for definition of parameters
        double dep = x[0];
        double age = p[0];
        double rho = p[1];
        double delta = p[2];

        // Parameters in Cherenkov yield
        double k_out = 2 * Pi() * Utility::fine_struct / rho * (1 / lambda_min - 1 / lambda_max);
        double k_1 = k_out * 2 * delta;
        double k_2 = k_out * Sq(Utility::mass_e);

        // Parameters in energy distribution
        double a1 = fe_a11 - fe_a12 * age;
        double a2 = fe_a21 - fe_a22 * age;
        double a0 = fe_k0 * Exp(fe_k1 * age + fe_k2 * Sq(age));

        // Product of the yield and energy distribution, integrating over lnE
        return a0 * Exp(dep) / ((a1 + Exp(dep)) * Power(a2 + Exp(dep), age)) * (k_1 - k_2 * Exp(-2.0 * dep));
    }

    void Simulator::ViewFluorescencePhotons(Shower shower, PhotonCount& photon_count)
    {
        int number_detected = NumberFluorescencePhotons(shower);
        for (int i = 0; i < number_detected / fluor_thin; i++)
        {
            TVector3 lens_impact = rotate_to_world * RandomStopImpact();
            Ray photon = Ray(shower.Position(), lens_impact - shower.Position(), JitteredTime(shower));
            photon.PropagateToPoint(lens_impact);
            SimulateOptics(photon, photon_count, fluor_thin);
        }
    }

    void Simulator::ViewCherenkovPhotons(Shower shower, Plane ground_plane, PhotonCount& photon_count)
    {
        int number_detected = NumberCherenkovPhotons(shower);
        for (int i = 0; i < number_detected / ckv_thin; i++)
        {
            Ray photon = GenerateCherenkovPhoton(shower);
            photon.PropagateToPlane(ground_plane);
            TVector3 stop_impact = rotate_to_world * RandomStopImpact();
            photon.PropagateToPoint(stop_impact);
            SimulateOptics(photon, photon_count, ckv_thin);
        }
    }

    int Simulator::NumberFluorescencePhotons(Shower shower)
    {
        // Find the yield via formula from Kakimoto
        double rho = shower.LocalRho();
        double term_1 = fluor_a1 / (1 + fluor_b1 * rho * Sqrt(atm_temp));
        double term_2 = fluor_a2 / (1 + fluor_b2 * rho * Sqrt(atm_temp));
        double yield = IonizationLossRate(shower) / dep_1_4 * (term_1 + term_2);

        // Find the number produced and the fraction captured
        double total = yield * shower.GaisserHillas() * depth_step;
        double fraction = SphereFraction(shower.Position()) * DetectorEfficiency();
        return Utility::RandomRound(total * fraction, rng);
    }

    int Simulator::NumberCherenkovPhotons(Shower shower)
    {
        // Perform the integration, ignoring portion above the shower energy
        ckv_integrator.SetParameter("age", shower.Age());
        ckv_integrator.SetParameter("rho", shower.LocalRho());
        ckv_integrator.SetParameter("delta", shower.LocalDelta());
        double yield = ckv_integrator.Integral(Log(shower.EThresh()), Log(shower.EnergyMeV()));

        // Find the number produced, and use a Lambertian reflectance model to find the fraction captured
        double total = yield * shower.GaisserHillas() * depth_step;
        TVector3 ground_impact = shower.PlaneImpact(ground_plane);
        double cos_theta = Abs(Cos(ground_impact.Angle(ground_plane.Normal())));
        double fraction = 4 * SphereFraction(ground_impact) * cos_theta * DetectorEfficiency();
        return Utility::RandomRound(total * fraction, rng);
    }

    void Simulator::SimulateOptics(Ray photon, PhotonCount& photon_count, double thinning)
    {
        // Refract across the corrector plate
        photon.Transform(rotate_to_world.Inverse());
        if (!DeflectFromLens(photon)) return;

        // Check for collision with the back of the pmt array
        TVector3 camera_impact;
        if (CameraImpactPoint(photon, camera_impact)) return;

        // Reflect the photon from the mirror
        TVector3 reflect_point;
        if (!MirrorImpactPoint(photon, reflect_point)) return;
        photon.PropagateToPoint(reflect_point);
        photon.Reflect(MirrorNormal(reflect_point));

        // Detection by the cluster
        if (!CameraImpactPoint(photon, camera_impact)) return;
        photon.PropagateToPoint(camera_impact);
        photon_count.AddPhoton(photon.Time(), camera_impact, thinning);
    }

    TVector3 Simulator::RandomStopImpact()
    {
        double r_rand = Utility::RandLinear(rng, stop_diameter / 2.0);
        double phi_rand = rng.Uniform(TwoPi());
        return TVector3(r_rand * Cos(phi_rand), r_rand * Sin(phi_rand), 0);
    }

    bool Simulator::DeflectFromLens(Ray& photon)
    {
        // Cut out the inner portion of the corrector
        double x = photon.Position().X();
        double y = photon.Position().Y();
        double photon_axis_dist = Sqrt(Sq(x) + Sq(y));
        if (photon_axis_dist < stop_diameter / (2.0 * Sqrt(2))) return photon.Direction().Z() < 0;

        // Find the normal to the corrector
        double c_lead = 4.0 * (refrac_lens - 1) * Power(mirror_radius, 3);
        double z_norm = c_lead / (4.0 * (Sq(x) + Sq(y)) - stop_diameter);
        TVector3 norm = TVector3(-x, -y, z_norm).Unit();

        // Perform the refraction
        bool success = photon.Refract(norm, 1, refrac_lens);
        return success && photon.Refract(TVector3(0, 0, 1), refrac_lens, 1);
    }

    bool Simulator::MirrorImpactPoint(Ray ray, TVector3& point)
    {
        NegSphereImpact(ray, point, mirror_radius);
        return Utility::WithinXYDisk(point, mirror_size / 2.0) && point.Z() < 0.0;
    }

    TVector3 Simulator::MirrorNormal(TVector3 point)
    {
        return -point.Unit();
    }

    bool Simulator::CameraImpactPoint(Ray ray, TVector3& point)
    {
        NegSphereImpact(ray, point, mirror_radius / 2.0);
        return Utility::WithinXYDisk(point, cluster_diameter / 2.0) && point.Z() < 0;
    }

    double Simulator::IonizationLossRate(Shower shower)
    {
        double age = shower.Age();
        return ion_c1 / Power(ion_c2 + age, ion_c3) + ion_c4 + ion_c5 * age;
    }

    double Simulator::SphereFraction(TVector3 view_point)
    {
        TVector3 detector_axis = rotate_to_world * TVector3(0, 0, 1);
        double cosine = Cos(detector_axis.Angle(view_point));
        cosine = cosine < 0.0 ? 0.0 : cosine;
        double area_fraction = Sq(stop_diameter / 2.0) / (4.0 * view_point.Mag2());
        return area_fraction * cosine;
    }

    double Simulator::DetectorEfficiency()
    {
        return quantum_eff * mirror_reflect * filter_transmit;
    }

    Ray Simulator::GenerateCherenkovPhoton(Shower shower)
    {
        TVector3 direction = shower.Direction();
        TVector3 rotation_axis = Utility::RandNormal(shower.Velocity().Unit(), rng);
        direction.Rotate(rng.Exp(ThetaC(shower)), rotation_axis);
        return Ray(shower.Position(), direction, JitteredTime(shower));
    }

    double Simulator::ThetaC(Shower shower)
    {
        return ckv_k1 * Power(shower.EThresh(), ckv_k2);
    }

    double Simulator::JitteredTime(Shower shower)
    {
        double step_time = depth_step / shower.LocalRho() / Utility::c_cent;
        return shower.Time() + rng.Uniform(-0.5 * step_time, 0.5 * step_time);
    }

    bool Simulator::NegSphereImpact(Ray ray, TVector3& point, double radius)
    {
        // Find roots of the constraining polynomial
        double a = ray.Velocity().Mag2();
        double b = 2 * ray.Position().Dot(ray.Velocity());
        double c = ray.Position().Mag2() - Sq(radius);
        ROOT::Math::Polynomial poly = ROOT::Math::Polynomial(a, b, c);
        std::vector<double> roots = poly.FindRealRoots();

        // If there are no real roots there is no impact, otherwise find the rearmost impact
        if (roots.size() == 0)
        {
            point = TVector3();
            return false;
        }
        else if (roots.size() == 1)
        {
            point = ray.Position() + roots[0] * ray.Velocity();
            return true;
        }
        else
        {
            TVector3 point1 = ray.Position() + roots[0] * ray.Velocity();
            TVector3 point2 = ray.Position() + roots[1] * ray.Velocity();
            point = point1.Z() < 0 ? point1 : point2;
            return true;
        }
    }

    double Simulator::MinTime(Shower shower) {
        double time = shower.Time();
        time += shower.Position().Mag() / Utility::c_cent;
        return time;
    }

    double Simulator::MaxTime(Shower shower) {
        double time = shower.Time();
        time += shower.TimeToPlane(ground_plane);
        time += shower.PlaneImpact(ground_plane).Mag() * back_toler / Utility::c_cent;
        return time;
    }
}

// Simulator.cpp
//
// Author: Matthew Dutson
//
// Implementation of Simulator.h

#include <boost/property_tree/ptree.hpp>
#include <TMath.h>

#include "Simulator.h"

using namespace TMath;

using std::string;
using boost::property_tree::ptree;

namespace cherenkov_simulator
{
    Simulator::Simulator(const ptree& config)
    {
        // Behavior of the simulation
        depth_step = config.get<double>("simulation.depth_step");
        fluor_thin = config.get<int>("simulation.fluor_thin");
        ckv_thin = config.get<int>("simulation.ckv_thin");
        back_toler = config.get<double>("simulation.back_toler");

        // Surrounding setup and orientation
        TVector3 ground_normal = Utility::ToVector(config.get<string>("surroundings.ground_normal"));
        TVector3 ground_point = Utility::ToVector(config.get<string>("surroundings.ground_point"));
        ground_plane = Plane(ground_normal, ground_point);
        rotate_to_world = Utility::MakeRotation(config.get<double>("surroundings.elevation_angle"));

        // Properties of the detector
        mirror_radius = config.get<double>("detector.mirror_radius");
        stop_diameter = mirror_radius / (2.0 * config.get<double>("detector.f_number"));
        auto field_of_view = config.get<double>("detector.field_of_view");
        mirror_size = stop_diameter + 2.0 * mirror_radius * Tan(field_of_view / 2.0);
        cluster_diameter = mirror_radius * Sin(field_of_view / 2.0);

        // PhotonCount parameters
        count_params.bin_size = config.get<double>("simulation.time_bin");
        count_params.max_bytes = config.get<size_t>("simulation.max_bytes");
        count_params.n_pixels = config.get<size_t>("detector.n_pmt_across");
        count_params.linear_size = cluster_diameter / count_params.n_pixels;
        count_params.angular_size = count_params.linear_size / (mirror_radius / 2.0);

        // Miscellaneous member variables
        ckv_integrator = TF1("ckv_integrator", ckv_func, 0.0, Infinity(), 3);
        ckv_integrator.SetParNames("age", "rho", "delta");
    }

    PhotonCount Simulator::SimulateShower(Shower shower) const
    {
        // A lower bound on the time the first photon will reach the detector
        PhotonCount photon_count = PhotonCount(count_params, MinTime(shower), MaxTime(shower));
        TF1 integrator = TF1(ckv_integrator);

        // Step the shower along its path
        while (shower.TimeToPlane(ground_plane) > 0)
        {
            shower.IncrementDepth(depth_step);
            ViewFluorescencePhotons(shower, photon_count);
            ViewCherenkovPhotons(shower, ground_plane, photon_count, integrator);
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

    void Simulator::ViewFluorescencePhotons(Shower shower, PhotonCount& photon_count) const
    {
        int n_loops = NumberFluorescenceLoops(shower);
        for (int i = 0; i < n_loops / fluor_thin; i++)
        {
            TVector3 lens_impact = rotate_to_world * RandomStopImpact();
            Ray photon = JitteredRay(shower, lens_impact - shower.Position());
            photon.PropagateToPoint(lens_impact);
            SimulateOptics(photon, photon_count, fluor_thin);
        }
    }

    void Simulator::ViewCherenkovPhotons(Shower shower, Plane ground_plane, PhotonCount& photon_count, TF1 integrator) const
    {
        int n_loops = NumberCherenkovLoops(shower, integrator);
        for (int i = 0; i < n_loops; i++)
        {
            Ray photon = GenerateCherenkovPhoton(shower);
            photon.PropagateToPlane(ground_plane);
            TVector3 stop_impact = rotate_to_world * RandomStopImpact();
            photon.PropagateToPoint(stop_impact);
            SimulateOptics(photon, photon_count, ckv_thin);
        }
    }

    int Simulator::NumberFluorescenceLoops(Shower shower) const
    {
        // Find the yield via formula from Kakimoto
        double rho = shower.LocalRho();
        double term_1 = fluor_a1 / (1 + fluor_b1 * rho * Sqrt(atm_temp));
        double term_2 = fluor_a2 / (1 + fluor_b2 * rho * Sqrt(atm_temp));
        double yield = IonizationLossRate(shower) / dep_1_4 * (term_1 + term_2);

        // Find the number produced and the fraction captured
        double total = yield * shower.GaisserHillas() * depth_step;
        double fraction = SphereFraction(shower.Position()) * DetectorEfficiency();
        return Utility::RandomRound(total * fraction / (double) fluor_thin);
    }

    int Simulator::NumberCherenkovLoops(Shower shower, TF1 integrator) const
    {
        // Perform the integration, ignoring portion above the shower energy
        integrator.SetParameter("age", shower.Age());
        integrator.SetParameter("rho", shower.LocalRho());
        integrator.SetParameter("delta", shower.LocalDelta());
        double yield = integrator.Integral(Log(shower.EThresh()), Log(shower.EnergyMeV()));

        // Find the number produced, and use a Lambertian reflectance model to find the fraction captured
        double total = yield * shower.GaisserHillas() * depth_step;
        TVector3 ground_impact = shower.PlaneImpact(ground_plane);
        double cos_theta = Abs(Cos(ground_impact.Angle(ground_plane.Normal())));
        double fraction = 4.0 * SphereFraction(ground_impact) * cos_theta * DetectorEfficiency();
        return Utility::RandomRound(total * fraction / (double) ckv_thin);
    }

    void Simulator::SimulateOptics(Ray photon, PhotonCount& photon_count, int thinning) const
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

    TVector3 Simulator::RandomStopImpact() const
    {
        double r_rand = Utility::RandLinear(0.0, stop_diameter / 2.0);
        double phi_rand = gRandom->Uniform(TwoPi());
        return TVector3(r_rand * Cos(phi_rand), r_rand * Sin(phi_rand), 0);
    }

    bool Simulator::DeflectFromLens(Ray& photon) const
    {
        // Cut out the inner portion of the corrector
        double x = photon.Position().X();
        double y = photon.Position().Y();
        double photon_axis_dist = Sqrt(Sq(x) + Sq(y));
        if (photon_axis_dist < stop_diameter / (2.0 * Sqrt(2))) return photon.Direction().Z() < 0;

        // Find the normal to the corrector
        double z_norm = (refrac_lens - 1) * Power(mirror_radius, 3) / (Sq(x) + Sq(y));
        TVector3 norm = TVector3(-x, -y, z_norm).Unit();

        // Perform the refraction
        bool success = photon.Refract(norm, 1, refrac_lens);
        return success && photon.Refract(TVector3(0, 0, 1), refrac_lens, 1);
    }

    bool Simulator::MirrorImpactPoint(Ray ray, TVector3& point) const
    {
        NegSphereImpact(std::move(ray), point, mirror_radius);
        return Utility::WithinXYDisk(point, mirror_size / 2.0) && point.Z() < 0.0;
    }

    TVector3 Simulator::MirrorNormal(TVector3 point) const
    {
        return -point.Unit();
    }

    bool Simulator::CameraImpactPoint(Ray ray, TVector3& point) const
    {
        NegSphereImpact(std::move(ray), point, mirror_radius / 2.0);
        return Utility::WithinXYDisk(point, cluster_diameter / 2.0) && point.Z() < 0;
    }

    double Simulator::IonizationLossRate(Shower shower) const
    {
        double age = shower.Age();
        return ion_c1 / Power(ion_c2 + age, ion_c3) + ion_c4 + ion_c5 * age;
    }

    double Simulator::SphereFraction(TVector3 view_point) const
    {
        TVector3 detector_axis = rotate_to_world * TVector3(0, 0, 1);
        double cosine = Cos(detector_axis.Angle(view_point));
        cosine = cosine < 0.0 ? 0.0 : cosine;
        double area_fraction = Sq(stop_diameter / 2.0) / (4.0 * view_point.Mag2());
        return area_fraction * cosine;
    }

    double Simulator::DetectorEfficiency() const
    {
        return quantum_eff * mirror_reflect * filter_transmit;
    }

    Ray Simulator::GenerateCherenkovPhoton(Shower shower) const
    {
        TVector3 direction = shower.Direction();
        TVector3 rotation_axis = Utility::RandNormal(shower.Velocity().Unit());
        direction.Rotate(gRandom->Exp(ThetaC(shower)), rotation_axis);
        return JitteredRay(shower, direction);
    }

    double Simulator::ThetaC(Shower shower) const
    {
        return ckv_k1 * Power(shower.EThresh(), ckv_k2);
    }

    Ray Simulator::JitteredRay(Shower shower, TVector3 direction) const
    {
        double step_time = depth_step / shower.LocalRho() / Utility::c_cent;
        double offset = gRandom->Uniform(-0.5 * step_time, 0.5 * step_time);
        double time = shower.Time() + offset;
        TVector3 position = shower.Position() + shower.Velocity() * offset;
        return Ray(position, direction, time);
    }

    bool Simulator::NegSphereImpact(Ray ray, TVector3& point, double radius)
    {
        // Find roots of the constraining polynomial
        double a = ray.Velocity().Mag2();
        double b = 2 * ray.Position().Dot(ray.Velocity());
        double c = ray.Position().Mag2() - Sq(radius);

        // If there are no real roots there is no impact, otherwise find the rearmost impact
        double b4ac = Sq(b) - 4 * a * c;
        if (b4ac < 0)
        {
            point = TVector3();
            return false;
        }
        else
        {
            double root1 = (-b + Sqrt(b4ac)) / (2 * a);
            double root2 = (-b - Sqrt(b4ac)) / (2 * a);
            TVector3 point1 = ray.Position() + root1 * ray.Velocity();
            TVector3 point2 = ray.Position() + root2 * ray.Velocity();
            point = point1.Z() < point2.Z() ? point1 : point2;
            return true;
        }
    }

    double Simulator::MinTime(Shower shower) const
    {
        double time = shower.Time();
        time += shower.Position().Mag() / Utility::c_cent;
        return time;
    }

    double Simulator::MaxTime(Shower shower) const
    {
        double time = shower.Time();
        time += shower.TimeToPlane(ground_plane);
        time += shower.PlaneImpact(ground_plane).Mag() * back_toler / Utility::c_cent;
        return time;
    }
}

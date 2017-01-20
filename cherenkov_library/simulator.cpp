// shower_simulator.cpp
// cherenkov_library
//
// Created by Matthew Dutson on 9/8/16.
//
// Contains the implementation of methods in simulator.h

#include <boost/property_tree/ptree.hpp>
#include <Math/Polynomial.h>

#include "simulator.h"
#include "utility.h"

using namespace TMath;

using std::string;
using boost::property_tree::ptree;

namespace cherenkov_library
{
    void Simulator::ParseFile(ptree config)
    {
        // Parameters related to the behavior of the simulation
        depth_step = config.get<double>("depth_step");
        time_bin = config.get<double>("time_bin");

        // Parameters relating to the position and orientation of the detector relative to its surroundings. The frames
        // share the same x-axis.
        ground_plane = Plane(ToVector(config.get<string>("ground_normal")),
                             ToVector(config.get<string>("ground_point")));
        rotate_to_world = TRotation();
        rotate_to_world.RotateX(-PiOver2() + config.get<double>("elevation_angle"));

        // The constant temperature of the atmosphere
        atmosphere_temp = config.get<double>("atmosphere_temp");

        // Parameters defining properties of the detector optics
        refrac_lens = config.get<double>("refrac_lens");
        mirror_radius = config.get<double>("mirror_radius");
        stop_size = config.get<double>("stop_size");
        mirror_size = config.get<double>("mirror_size");
        cluster_size = config.get<double>("cluster_size");
        n_pmt_across = config.get<int>("n_pmt_across");

        // Parameters in the GH profile
        gh_lambda = config.get<double>("gh_lambda");

        // Parameters used when calculating the fluorescence yield
        fluor_a1 = config.get<double>("fluor_a1");
        fluor_a2 = config.get<double>("fluor_a2");
        fluor_b1 = config.get<double>("fluor_b1");
        fluor_b2 = config.get<double>("fluor_b2");
        dep_1_4 = config.get<double>("dep_1_4");

        // Parameters used when calculating the effective ionization loss rate
        ion_c1 = config.get<double>("ion_c1");
        ion_c2 = config.get<double>("ion_c2");
        ion_c3 = config.get<double>("ion_c3");
        ion_c4 = config.get<double>("ion_c4");
        ion_c5 = config.get<double>("ion_c5");

        // Parameters used when calculating theta_c in the Cherenkov angular distribution
        ckv_k1 = config.get<double>("ckv_k1");
        ckv_k2 = config.get<double>("ckv_k2");

        // Parameters used when calculating the Cherenkov yield
        lambda_min = config.get<double>("lambda_min");
        lambda_max = config.get<double>("lambda_max");

        // Parameters in the electrohn energy spectrum
        fe_a11 = config.get<double>("fe_a11");
        fe_a12 = config.get<double>("fe_a12");
        fe_a21 = config.get<double>("fe_a21");
        fe_a22 = config.get<double>("fe_a22");
        fe_k0 = config.get<double>("fe_k0");
        fe_k1 = config.get<double>("fe_k1");
        fe_k2 = config.get<double>("fe_k2");

        // Physics constants
        mass_e = config.get<double>("mass_e");
        fine_struct = config.get<double>("fine_struct");

        // Parameters defining the amount of night sky background noise
        sky_noise = config.get<double>("sky_noise");
        ground_noise = config.get<double>("ground_noise");

        // A general-purpose random number generator
        rng = TRandom3();
    }

    PhotonCount Simulator::SimulateShower(Shower shower)
    {
        // Arc length = angle * radius
        double pmt_linear_size = cluster_size / n_pmt_across;
        double pmt_angular_size = pmt_linear_size / (mirror_radius / 2.0);
        PhotonCount photon_count = PhotonCount(n_pmt_across, shower.Time(), time_bin, pmt_angular_size,
                                               pmt_linear_size);

        // Step the shower through its path.
        while (shower.TimeToPlane(ground_plane) > 0)
        {
            shower.IncrementDepth(depth_step);
            ViewFluorescencePhotons(shower, &photon_count);
            ViewCherenkovPhotons(shower, ground_plane, &photon_count);
        }

        AddNoise(&photon_count);
        return photon_count;
    }

    void Simulator::ViewFluorescencePhotons(Shower shower, PhotonCount* photon_count)
    {
        int number_detected = NumberFluorescencePhotons(shower);

        for (int i = 0; i < number_detected; i++)
        {
            TVector3 lens_impact = RandomStopImpact();
            Ray photon = Ray(shower.Position(), lens_impact - shower.Position(), shower.Time());
            photon.PropagateToPoint(lens_impact);

            SimulateOptics(photon, photon_count);
        }
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

    TVector3 Simulator::RandomStopImpact()
    {
        double phi_rand = rng.Uniform(TwoPi());

        // The probability that a random point on a circle will lie at radius r is proportional to r (area of slice is
        // 2*pi*r*dr).
        TF1 radius_distribution = TF1("stop_impact", "x", 0, stop_size / 2.0);
        double r_rand = radius_distribution.GetRandom();
        return TVector3(r_rand * Cos(phi_rand), r_rand * Sin(phi_rand), 0);
    }

    void Simulator::DeflectFromLens(Ray* photon)
    {
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

        // TODO: Here we are assuming an index of refraction for the air of one. Should we change this?
        photon->Refract(norm, 1, refrac_lens);

        photon->Refract(TVector3(0, 0, 1), refrac_lens, 1);

    }

    bool Simulator::MirrorImpactPoint(Ray ray, TVector3* point)
    {
        // Find the point with the smallest z-component where the ray intersects with the mirror sphere
        BackOriginSphereImpact(ray, point, mirror_radius);

        // Ensure that the ray actually hits the mirror
        return WithinXYDisk(*point, mirror_size / 2.0) && point->Z() < 0.0;
    }

    TVector3 Simulator::MirrorNormal(TVector3 point)
    {
        // For a spherical mirror, the normal vector will always point straight back to the center of curvature (the
        // origin in this case).
        return -point;
    }

    bool Simulator::CameraImpactPoint(Ray ray, TVector3* point)
    {
        BackOriginSphereImpact(ray, point, mirror_radius);
        return WithinXYDisk(*point, cluster_size / 2.0) && point->Z() < 0;
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
        }
        else
        {
            if (ray.Velocity().Z() < 0)
            {
                time = Max(roots[0], roots[1]);
            }
            else
            {
                time = Min(roots[0], roots[1]);
            }
        }
        *point = ray.Position() + time * ray.Velocity();
        return true;
    }

    int Simulator::NumberFluorescencePhotons(Shower shower)
    {
        int n_charged = GaiserHilles(shower);
        double alpha_eff = IonizationLossRate(shower);

        double rho = shower.LocalRho();
        double term_1 = fluor_a1 / (1 + fluor_b1 * rho * Sqrt(atmosphere_temp));
        double term_2 = fluor_a2 / (1 + fluor_b2 * rho * Sqrt(atmosphere_temp));

        // This yield matches the form of Stratton 4.2 (from Kakimoto). The energy deposit rate for a single photon is
        // alpha_eff, so the deposit rate for all photons is alpha_eff * N.
        double yield = alpha_eff / dep_1_4 * rho * (term_1 + term_2);
        double total_produced = yield * n_charged * depth_step;

        // Find the fraction captured by the camera.
        return total_produced * SphereFraction(shower.Position());
    }

    int Simulator::NumberCherenkovPhotons(Shower shower)
    {
        // 2 * pi * alpha * (1 / lambda1 - 1 / lambda2) / rho
        double rho = shower.LocalRho();
        double k_out = 2 * Pi() * fine_struct / rho * (1 / lambda_min - 1 / lambda_max);
        double k_1 = k_out * 2 * shower.LocalDelta();

        // We don't need to multiply by c^2 here because our electron mass should already be in units of MeV/c^2.
        double k_2 = k_out * Sq(mass_e);

        // Parameters in the electron energy distribution
        double age = shower.Age();
        double a1 = fe_a11 - fe_a12 * age;
        double a2 = fe_a21 - fe_a22 * age;
        double a0 = fe_k0 * Exp(fe_k1 * age + fe_k2 * Sq(age));

        // See notes for details. Integrate over lnE from lnE_thresh to infinity. E terms turn to e^E because the
        // variable of integration is lnE.
        std::stringstream func_string;
        func_string << "(" << a0 << "* e^x / ((" << a1 << "+e^x)(" << a2 << "+e^x)^" << age << "))(" << k_1 << "-"
                    << k_2 << "/e^x)";
        double e_thresh = EThresh(shower);
        TF1 func = TF1("integrand", func_string.str().c_str(), e_thresh, Infinity());
        double integral = func.Integral(e_thresh, Infinity());

        // Determine the number captured. Multiply the PhotonFraction by two because we're dealing with a half sphere,
        // not a full sphere.
        return GaiserHilles(shower) * integral * depth_step * SphereFraction(shower.PlaneImpact(ground_plane)) * 2;
    }

    double Simulator::GaiserHilles(Shower shower)
    {
        // Use the Gaisser-Hilles profile to find the number of electrons.
        double x = shower.X();
        double n_max = shower.NMax();
        double x_max = shower.XMax();
        double x_0 = shower.X0();
        return n_max * Power((x - x_0) / (x_max - x_0), (x_max - x_0) / gh_lambda) * Exp((x_max - x) / gh_lambda);
    }

    double Simulator::IonizationLossRate(Shower shower)
    {
        double s = shower.Age();
        return ion_c1 / Power(ion_c2 + s, ion_c3) + ion_c4 + ion_c5 * s;
    }

    double Simulator::EThresh(Shower shower)
    {
        double delta = shower.LocalDelta();

        // We don't need to multiply by c^2 here because our electron mass should already be in units of MeV/c^2.
        return mass_e / Sqrt(2 * delta);
    }

    double Simulator::SphereFraction(TVector3 view_point)
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
        std::string formula = "e^(-x/" + std::to_string(ThetaC(shower)) + ")/sin(x)";
        TF1 angular_distribution = TF1("distro", formula.c_str(), 0, Pi());
        TVector3 direction = shower.Velocity().Unit();
        direction.Rotate(angular_distribution.GetRandom(), rotation_axis);
        return Ray(shower.Position(), direction, shower.Time());
    }

    double Simulator::ThetaC(Shower shower)
    {
        return ckv_k1 * Power(EThresh(shower), ckv_k2);
    }

    void Simulator::AddNoise(PhotonCount* photon_count)
    {
        SignalIterator iter = photon_count->Iterator();
        while (iter.Next())
        {
            TVector3 direction = photon_count->Direction(iter);
            Ray outward_ray = Ray(direction, TVector3(0, 0, 0), 0);

            // If the pixel is looking at the ground, use the noise rate for the ground. Otherwise, use the sky rate.
            double noise_rate;
            if (outward_ray.TimeToPlane(ground_plane) > 0)
            {
                noise_rate = ground_noise;
            }
            else
            {
                noise_rate = sky_noise;
            }
            photon_count->AddNoise(noise_rate, iter, rng);
        }
    }
}

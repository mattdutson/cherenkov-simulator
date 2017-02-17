// reconstruction.cpp
// cherenkov_library
//
// Created by Matthew Dutson on 9/8/16.
//
//

#include "reconstructor.h"
#include "data_containers.h"
#include "TMatrixD.h"
#include "TMatrixDSymEigen.h"
#include "TVectorD.h"
#include "utility.h"
#include <string>
#include "TMath.h"
#include <Math/Functor.h>
#include <Math/GSLMinimizer.h>
#include <TGraphErrors.h>
#include <TF1.h>

using boost::property_tree::ptree;
using std::string;
using std::vector;
using namespace TMath;

namespace cherenkov_library
{
    Reconstructor::Reconstructor(ptree config)
    {
        // Construct the ground plane in the world frame.
        ground_plane = Plane(ToVector(config.get<string>("ground_normal")),
                             ToVector(config.get<string>("ground_point")));

        // The rotation from detector frame to world frame. The frames share the same x-axis.
        rotate_to_world = TRotation();
        rotate_to_world.RotateX(-PiOver2() + config.get<double>("elevation_angle"));

        // Parameters defining the amount of night sky background noise
        double stop_diameter = config.get<double>("mirror_radius") / (2.0 * config.get<double>("f_number"));
        sky_noise = config.get<double>("sky_noise") * Pi() * Sq(stop_diameter / 2.0);
        ground_noise = config.get<double>("ground_noise") * Pi() * Sq(stop_diameter / 2.0);

        // Parameters used when applying triggering logic
        trigger_thresh = config.get<double>("trigger_thresh");
        hold_thresh = config.get<double>("hold_thresh");
        trigger_clust = config.get<int>("trigger_clust");
    }

    Plane Reconstructor::FitSDPlane(PhotonCount data)
    {
        // Compute the matrix specified in Stratton 3.4. Maybe just iterate over matrix elements whose symmetric
        // counterparts haven't been found yet.
        SignalIterator iter = data.Iterator();
        TMatrixDSym matrix(3);
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                double mat_element = 0;
                iter.Reset();
                while (iter.Next())
                {
                    TVector3 direction = data.Direction(iter);
                    int pmt_sum = data.SumBins(iter);
                    mat_element += direction[j] * direction[k] * pmt_sum;
                }
                matrix[j][k] = mat_element;
            }
        }

        // Find the eigenvector corresponding to the minimum eigenvalue.
        TMatrixDSymEigen eigen = TMatrixDSymEigen(matrix);
        TVectorD eigen_val = eigen.GetEigenValues();
        TMatrixD eigen_vec = eigen.GetEigenVectors();
        double min_val = eigen_val[0];
        int min_index = 0;
        for (int i = 0; i < eigen_val.GetNoElements(); i++)
        {
            if (eigen_val[i] < min_val)
            {
                min_val = eigen_val[i];
                min_index = i;
            }
        }

        // Each eigenvector is stored as a column in the matrix. Matrices are indexed as [row][column].
        TVector3 result = TVector3(eigen_vec[0][min_index], eigen_vec[1][min_index], eigen_vec[2][min_index]);
        return Plane(rotate_to_world * result, TVector3());
    }

    TGraph
    Reconstructor::MonocularFit(PhotonCount data, Plane sd_plane, double* t_0, double* impact_param, double* angle)
    {
        TGraphErrors data_graph = GetFitGraph(data, sd_plane);

        // The functional form of the time profile.
        std::stringstream func_string = std::stringstream();
        func_string << "[0] + [1] / (" << CentC() << ") * tan((" << Pi() << " - [2] - x) / 2)";
        TF1 func = TF1("profile_fit", func_string.str().c_str(), -Pi(), Pi());

        // Set names and initial guesses for parameters.
        func.SetParNames("t_0", "r_p", "psi");
        func.SetParameters(0.0, 10e6, PiOver2());

        // Perform the fit.
        data_graph.Fit("profile_fit");
        TF1* result = data_graph.GetFunction("profile_fit");
        *t_0 = result->GetParameter("t_0");
        *impact_param = result->GetParameter("r_p");
        *angle = result->GetParameter("psi");

        // Return a graph containing the data points with error bars.
        return data_graph;
    }

    TGraph
    Reconstructor::HybridFit(PhotonCount data, TVector3 impact_point, Plane sd_plane, double* t_0, double* impact_param,
                             double* angle)
    {
        TGraphErrors data_graph = GetFitGraph(data, sd_plane);

        // Find the angle defining the relationship between the impact angle and impact parameter.
        double impact_distance = impact_point.Mag();
        TVector3 horizontal = TVector3(0, 1, 0);
        TVector3 projection = (horizontal - horizontal.Dot(sd_plane.Normal()) * sd_plane.Normal()).Unit();
        double theta = impact_point.Angle(projection);
        if (Above(horizontal, impact_point))
        {
            theta = Abs(theta);
        }
        else
        {
            theta = -Abs(theta);
        }

        // The functional form of the time profile.
        std::stringstream func_string = std::stringstream();
        func_string << "[0] + " << impact_distance << " * sin([1] - " << theta << ") / (" << CentC() << ") * tan(("
                    << Pi() << " - [1] - x) / 2)";
        TF1 func = TF1("profile_fit", func_string.str().c_str(), -Pi(), Pi());

        // Set names and initial guesses for parameters.
        func.SetParNames("t_0", "psi");
        func.SetParameters(0.0, PiOver2());

        // Perform the fit.
        data_graph.Fit("profile_fit");
        TF1* result = data_graph.GetFunction("profile_fit");
        *t_0 = result->GetParameter("t_0");
        *angle = result->GetParameter("psi");
        *impact_param = impact_distance * Sin(*angle);

        // Return a graph containing the data points with error bars.
        return data_graph;
    }

    bool Reconstructor::ApplyTriggering(PhotonCount* data)
    {
        // Initialize the structure to contain triggering values
        int size = data->Size();
        vector<vector<vector<bool>>> triggers = vector<vector<vector<bool>>>(size, vector<vector<bool>>(size));

        SignalIterator iter = data->Iterator();
        while (iter.Next())
        {
            // Determine whether a pixel is looking at the ground or sky.
            int noise_rate;
            Ray outward_ray = Ray(TVector3(), rotate_to_world * data->Direction(iter), 0);
            if (outward_ray.TimeToPlane(ground_plane) > 0)
            {
                noise_rate = ground_noise;
            }
            else
            {
                noise_rate = sky_noise;
            }

            // Remove any below-threshold noise in each pixel.
            data->ClearNoise(iter, noise_rate, hold_thresh);

            // Find any trigger times in the pixel.
            triggers[iter.X()][iter.Y()] = data->FindTriggers(iter, noise_rate, trigger_thresh);
        }

        // Determine which bins contain a sufficiently large cluster of triggered tubes
        vector<bool> good_frames = vector<bool>();
        bool triggered_found = false;
        for (int i = 0; i < data->NBins(); i++)
        {
            // Copy the true/false triggering values into a structure to be passed to LargestCluster.
            vector<vector<bool>> frame = vector<vector<bool>>(size, vector<bool>(size, false));
            for (int x = 0; x < triggers.size(); x++)
            {
                for (int y = 0; y < triggers[x].size(); y++)
                {
                    frame[x][y] = triggers[x][y][i];
                }
            }

            // Check whether there are any sufficiently large clusters of triggered tubes.
            if (LargestCluster(frame) > trigger_clust)
            {
                triggered_found = true;
                good_frames.push_back(true);
            }
            else
            {
                good_frames.push_back(false);
            }
        }

        // Zero all time bins in which there was no system trigger.
        data->EraseNonTriggered(good_frames);

        // Return whether any time bins contained a trigger.
        return triggered_found;
    }

    Shower Reconstructor::Reconstruct(PhotonCount data, bool try_ground, bool* triggered, bool* ground_used)
    {
        SubtractNoise(&data);
        *triggered = ApplyTriggering(&data);
        if (*triggered)
        {
            Plane sd_plane = FitSDPlane(data);
            double t_0, impact_param, angle;
            if (try_ground)
            {
                TVector3 impact;
                try_ground = FindGroundImpact(data, &impact);
                if (try_ground)
                {
                    HybridFit(data, impact, sd_plane, &t_0, &impact_param, &angle);
                }
            }
            if (!try_ground)
            {
                MonocularFit(data, sd_plane, &t_0, &impact_param, &angle);
            }

            // The direction in the shower-detector plane, with the y-axis defined to lie in the world's xy plane.
            // TODO: Check the signs on this transformation, as well as the direction of sd_plane.Normal x impact_direction
            TVector3 impact_direction = TVector3(Sin(angle), Cos(angle), 0);
            impact_direction.RotateUz(sd_plane.Normal());
            return Shower(Shower::Params(), sd_plane.Normal().Cross(impact_direction), impact_param * impact_direction,
                          t_0);
        }
        else
        {
            return Shower(Shower::Params(), TVector3(), TVector3());
        }
    }

    TGraphErrors Reconstructor::GetFitGraph(PhotonCount data, Plane sd_plane)
    {
        // Populate a TGraph with angle/time points from the data. Set errors based on the number of photons viewed.
        vector<double> angles = vector<double>();
        vector<double> times = vector<double>();
        vector<double> time_err = vector<double>();
        SignalIterator iter = data.Iterator();
        while (iter.Next())
        {
            TVector3 direction = rotate_to_world * data.Direction(iter);

            // Only consider pixels with a nonzero signal and those which point above the horizon.
            Ray outward_ray = Ray(TVector3(), direction, 0);
            if (outward_ray.TimeToPlane(ground_plane) < 0 && data.SumBins(iter) > 0)
            {
                // Project the direction vector onto the shower-detector plane and find the angle of the pixel.
                TVector3 projection = (direction - direction.Dot(sd_plane.Normal()) * sd_plane.Normal()).Unit();
                TVector3 horizontal = TVector3(0, 1, 0);

                // Make sure the angle has the correct sign.
                double chi = projection.Angle(horizontal);
                if (Above(horizontal, projection))
                {
                    chi = Abs(chi);
                }
                else
                {
                    chi = -Abs(chi);
                }

                // Add data points to the arrays.
                angles.push_back(chi);
                times.push_back(data.AverageTime(iter));
                time_err.push_back(1.0 / Sqrt((double) data.SumBins(iter)));
            }
        }

        // Make arrays to pass to the TGraph constructor.
        vector<double> angle_err = vector<double>(angles.size(), 0.0);
        return TGraphErrors(angles.size(), &angles[0], &times[0], &angle_err[0], &time_err[0]);
    }

    int Reconstructor::LargestCluster(std::vector<std::vector<bool>> not_counted)
    {
        int largest = 0;
        for (int i = 0; i < not_counted.size(); i++)
        {
            for (int j = 0; j < not_counted[i].size(); j++)
            {
                if (not_counted[i][j])
                {
                    int size = Visit(i, j, &not_counted);
                    if (size > largest)
                    {
                        largest = size;
                    }
                }
            }
        }
        return largest;
    }

    int Reconstructor::Visit(int i, int j, std::vector<std::vector<bool>>* not_counted)
    {
        if (i > not_counted->size() - 1 || i < 0)
        {
            return 0;
        }
        else if (j > not_counted->at(i).size() - 1 || j < 0)
        {
            return 0;
        }
        else if (!not_counted->at(i)[j])
        {
            return 0;
        }
        else
        {
            not_counted->at(i)[j] = false;
            int count = 1;
            count += Visit(i + 1, j, not_counted);
            count += Visit(i + 1, j + 1, not_counted);
            count += Visit(i - 1, j, not_counted);
            count += Visit(i - 1, j - 1, not_counted);
            count += Visit(i, j + 1, not_counted);
            count += Visit(i - 1, j + 1, not_counted);
            count += Visit(i, j - 1, not_counted);
            count += Visit(i + 1, j - 1, not_counted);
            return count;
        }
    }

    void Reconstructor::SubtractNoise(PhotonCount* data)
    {
        SignalIterator iter = data->Iterator();
        while (iter.Next())
        {
            Ray outward_ray = Ray(TVector3(), rotate_to_world * data->Direction(iter), 0);
            if (outward_ray.TimeToPlane(ground_plane) > 0)
            {
                data->SubtractNoise(ground_noise, iter);
            }
            else
            {
                data->SubtractNoise(sky_noise, iter);
            }
        }
    }

    bool Reconstructor::FindGroundImpact(PhotonCount data, TVector3* impact)
    {
        // Find the brightest pixel below the horizon.
        TVector3 impact_direction = TVector3();
        int highest_count = 0;
        SignalIterator iter = data.Iterator();
        while (iter.Next())
        {
            int sum = data.SumBins(iter);
            if (sum > highest_count)
            {
                highest_count = sum;
                impact_direction = rotate_to_world * data.Direction(iter);
            }
        }

        // Find the impact position by extending the direction to the ground.
        Ray outward_ray = Ray(TVector3(), impact_direction, 0);
        outward_ray.PropagateToPlane(ground_plane);
        *impact = outward_ray.Position();

        // Ensure that the number of photons seen in this brightest pixel exceeds some threshold.
        double poisson = data.RealNoiseRate(ground_noise) * data.NBins();
        double threshold = 3 * Sqrt(poisson);
        return highest_count > threshold;
    }
}

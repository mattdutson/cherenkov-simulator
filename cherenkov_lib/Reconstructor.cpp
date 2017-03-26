// Reconstructor.cpp
// cherenkov_lib
//
// Author: Matthew Dutson
//
// Implementation of Reconstructor.h

#include <TMatrixDSymEigen.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TF1.h>

#include "Reconstructor.h"
#include "Utility.h"

using namespace TMath;

using std::string;
using std::vector;
using boost::property_tree::ptree;

namespace cherenkov_simulator
{
    Reconstructor::Reconstructor(ptree config)
    {
        // Construct the ground plane in the world frame.
        TVector3 ground_normal = Utility::ToVector(config.get<string>("surroundings.ground_normal"));
        TVector3 ground_point = Utility::ToVector(config.get<string>("surroundings.ground_point"));
        ground_plane = Plane(ground_normal, ground_point);

        // The rotation from detector frame to world frame. The frames share the same x-axis.
        rotate_to_world = Utility::MakeRotation(config.get<double>("surroundings.elevation_angle"));

        // The amount of night sky background noise
        double stop_diameter = config.get<double>("detector.mirror_radius") / (2.0 * config.get<double>("detector.f_number"));
        sky_noise = Sq(stop_diameter / 2.0) * Pi() * global_sky_noise;
        ground_noise = Sq(stop_diameter / 2.0) * Pi() * global_ground_noise;

        // Parameters used when applying triggering logic
        trigger_thresh = config.get<double>("triggering.trigger_thresh");
        noise_thresh = config.get<double>("triggering.noise_thresh");
        trigger_clust = config.get<int>("triggering.trigger_clust");

        impact_buffer = config.get<double>("triggering.impact_buffer");
        rng = TRandom3();
    }

    Shower Reconstructor::Reconstruct(PhotonCount data, bool try_ground, bool* triggered, bool* ground_used)
    {
        *triggered = DetectorTriggered(GetTriggeringState(GetThresholdMatrices(data, trigger_thresh)));
        *ground_used = false;
        if (*triggered)
        {
            TRotation to_sd_plane = FitSDPlane(data);
            double t_0 = 0.0, impact_param = 0.0, angle = 0.0;
            MonocularFit(data, to_sd_plane, &t_0, &impact_param, &angle);
            Shower output = MakeShower(t_0, impact_param, angle, to_sd_plane);
            try_ground = PointWithinView(rotate_to_world.Inverse() * output.PlaneImpact(ground_plane), data);
            if (try_ground)
            {
                TVector3 impact;
                try_ground = FindGroundImpact(data, &impact);
                if (try_ground)
                {
                    try_ground = PointWithinView(impact, data);
                    if (try_ground)
                    {
                        HybridFit(data, impact, to_sd_plane, &t_0, &impact_param, &angle);
                        output = MakeShower(t_0, impact_param, angle, to_sd_plane);
                        *ground_used = true;
                    }
                }
            }
            return output;
        }
        else
        {
            return Shower(Shower::Params(), TVector3(), TVector3());
        }
    }

    void Reconstructor::AddNoise(PhotonCount* data)
    {
        PhotonCount::Iterator iter = data->GetIterator();
        while (iter.Next())
        {
            bool toward_ground = ground_plane.InFrontOf(rotate_to_world * data->Direction(&iter));
            data->AddNoise(toward_ground ? ground_noise : sky_noise, &iter, &rng);
        }
    }

    void Reconstructor::ClearNoise(PhotonCount* data)
    {
        SubtractAverageNoise(data);
        ThreeSigmaFilter(data);
        RecursiveSearch(data);
    }

    TGraphErrors
    Reconstructor::MonocularFit(PhotonCount data, TRotation to_sdp, double* t_0, double* r_p, double* psi)
    {
        TGraphErrors data_graph = GetFitGraph(data, to_sdp);

        // The functional form of the time profile.
        std::stringstream func_string = std::stringstream();
        func_string << "[0] + [1] / (" << Utility::c_cent << ") * tan((pi - [2] - x) / 2)";
        TF1 func = TF1("profile_fit", func_string.str().c_str(), -Pi(), Pi());

        // Set names and initial guesses for parameters.
        func.SetParNames("t_0", "r_p", "psi");
        func.SetParameters(0.0, 1e6, PiOver2());

        // Perform the fit.
        data_graph.Fit("profile_fit");
        TF1* result = data_graph.GetFunction("profile_fit");
        *t_0 = result->GetParameter("t_0");
        *r_p = result->GetParameter("r_p");
        *psi = result->GetParameter("psi");

        // Return a graph containing the data points with error bars.
        return data_graph;
    }

    TGraphErrors
    Reconstructor::HybridFit(PhotonCount data, TVector3 impact, TRotation to_sdp, double* t_0, double* r_p, double* psi)
    {
        TGraphErrors data_graph = GetFitGraph(data, to_sdp);

        // Find the angle of the impact direction with the shower-detector frame x-axis.
        double impact_distance = impact.Mag();
        double theta = (to_sdp * impact).Phi();

        // The functional form of the time profile.
        std::stringstream func_string = std::stringstream();
        func_string << "[0] +" << impact_distance << " * sin([1] + " << theta << ") / (" << Utility::c_cent << ") * tan((pi - [1] - x) / 2)";
        TF1 func = TF1("profile_fit", func_string.str().c_str(), -Pi(), Pi());

        // Set names and initial guesses for parameters.
        func.SetParNames("t_0", "psi");
        func.SetParameters(0.0, PiOver2());

        // Perform the fit.
        data_graph.Fit("profile_fit");
        TF1* result = data_graph.GetFunction("profile_fit");
        *t_0 = result->GetParameter("t_0");
        *psi = result->GetParameter("psi");
        *r_p = impact_distance * Sin(*psi);

        // Return a graph containing the data points with error bars.
        return data_graph;
    }

    TRotation Reconstructor::FitSDPlane(PhotonCount data)
    {
        // Compute the matrix specified in Stratton 3.4. Maybe just iterate over matrix elements whose symmetric
        // counterparts haven't been found yet.
        PhotonCount::Iterator iter = data.GetIterator();
        TMatrixDSym matrix(3);
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                double mat_element = 0;
                iter.Reset();
                while (iter.Next())
                {
                    TVector3 direction = data.Direction(&iter);
                    int pmt_sum = data.SumBins(&iter);
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
        TVector3 normal = TVector3(eigen_vec[0][min_index], eigen_vec[1][min_index], eigen_vec[2][min_index]);
        normal = rotate_to_world * normal;
        if (normal.X() < 0) normal = -normal;

        // Construct a rotation which takes points to the frame where the shower-detector plane is the xy plane.
        TVector3 new_x = (normal == TVector3(0, 0, 1)) ? TVector3(1, 0, 0) : TVector3(0, 0, 1).Cross(normal).Unit();
        TVector3 new_y = normal.Cross(new_x).Unit();
        return TRotation().RotateAxes(new_x, new_y, normal).Inverse();
    }

    bool Reconstructor::FindGroundImpact(PhotonCount data, TVector3* impact)
    {
        // Find the brightest pixel below the horizon.
        TVector3 impact_direction = TVector3();
        int highest_count = 0;
        PhotonCount::Iterator iter = data.GetIterator();
        while (iter.Next())
        {
            TVector3 direction = rotate_to_world * data.Direction(&iter);
            int sum = data.SumBins(&iter);
            if (sum > highest_count && ground_plane.InFrontOf(direction))
            {
                highest_count = sum;
                impact_direction = direction;
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

    TGraphErrors Reconstructor::GetFitGraph(PhotonCount data, TRotation to_sdp)
    {
        // Populate a TGraph with angle/time points from the data. Set errors based on the number of photons viewed.
        vector<double> angles = vector<double>();
        vector<double> times = vector<double>();
        vector<double> time_er = vector<double>();
        PhotonCount::Iterator iter = data.GetIterator();
        while (iter.Next())
        {
            // Don't rotate to the world because the rotation goes from the detector frame to the shower-detector frame.
            TVector3 direction = rotate_to_world * data.Direction(&iter);
            int bin_sum = data.SumBins(&iter);
            if (!ground_plane.InFrontOf(direction) && bin_sum > 0)
            {
                angles.push_back((to_sdp * direction).Phi());
                times.push_back(data.AverageTime(&iter));

                // Assume the error in each average time is the error in an individual time (bin size) over sqrt(n)
                time_er.push_back(data.BinSize() / Sqrt((double) bin_sum));
            }
        }

        // Make arrays to pass to the TGraph constructor.
        vector<double> angle_er = vector<double>(angles.size(), 0.0);
        TGraphErrors graph = TGraphErrors(angles.size(), &(angles[0]), &(times[0]), &(angle_er[0]), &(time_er[0]));
        graph.Sort();
        return graph;
    }

    void Reconstructor::SubtractAverageNoise(PhotonCount* data)
    {
        PhotonCount::Iterator iter = data->GetIterator();
        while (iter.Next())
        {
            bool toward_ground = ground_plane.InFrontOf(rotate_to_world * data->Direction(&iter));
            data->SubtractNoise(toward_ground ? ground_noise : sky_noise, &iter);
        }
    }

    void Reconstructor::ThreeSigmaFilter(PhotonCount* data)
    {
        PhotonCount::Iterator iter = data->GetIterator();
        while (iter.Next())
        {
            bool toward_ground = ground_plane.InFrontOf(rotate_to_world * data->Direction(&iter));
            data->ClearNoise(&iter, toward_ground ? ground_noise : sky_noise, noise_thresh);
        }
    }

    std::vector<bool> Reconstructor::GetTriggeringState(Bool3D trig_matrices)
    {
        vector<bool> good_frames = vector<bool>();
        for (int i = 0; i < trig_matrices[0][0].size(); i++) good_frames.push_back(FrameTriggered(i, &trig_matrices));
        return good_frames;
    }

    void Reconstructor::RecursiveSearch(PhotonCount* data)
    {
        // Find the pixels and times where triggers occured as well as pixels and times above the noise threshold
        Bool3D triggered = GetThresholdMatrices(*data, trigger_thresh);
        Bool3D three_sigma = GetThresholdMatrices(*data, noise_thresh);
        Bool3D good_pixels = data->GetFalseMatrix();
        vector<bool> trig_state = GetTriggeringState(GetThresholdMatrices(*data, trigger_thresh));

        // Start the recursive bleed from triggered pixels
        for (int i = 0; i < triggered.size(); i++)
        {
            for (int j = 0; j < triggered[i].size(); j++)
            {
                for (int t = 0; t < triggered[i][j].size(); t++)
                {
                    if (triggered[i][j][t] && trig_state[t])
                    {
                        BleedTrigger(i, j, t, &three_sigma, &good_pixels);
                    }
                }
            }
        }

        // Erase any pixels we determined to be non-valid
        data->Subset(good_pixels);
    }

    bool Reconstructor::DetectorTriggered(vector<bool> trig_state)
    {
        for (int i = 0; i < trig_state.size(); i++) if (trig_state[i]) return true;
        return false;
    }

    bool Reconstructor::FrameTriggered(int t, Bool3D* triggers)
    {
        return LargestCluster(t, triggers) > trigger_clust;
    }

    Reconstructor::Bool3D Reconstructor::GetThresholdMatrices(PhotonCount data, double sigma_mult, bool use_below_horiz)
    {
        Bool3D pass = data.GetFalseMatrix();
        PhotonCount::Iterator iter = data.GetIterator();
        while (iter.Next())
        {
            bool toward_ground = ground_plane.InFrontOf(rotate_to_world * data.Direction(&iter));
            if (toward_ground && !use_below_horiz) continue;
            pass[iter.X()][iter.Y()] = data.AboveThreshold(&iter, toward_ground ? ground_noise : sky_noise, sigma_mult);
        }
        return pass;
    }

    int Reconstructor::LargestCluster(int t, Bool3D* not_counted)
    {
        int largest = 0;
        for (int i = 0; i < not_counted->size(); i++)
        {
            for (int j = 0; j < not_counted->at(i).size(); j++)
            {
                if (not_counted->at(i)[j][t])
                {
                    int size = Visit(i, j, t, not_counted);
                    if (size > largest) largest = size;
                }
            }
        }
        return largest;
    }

    int Reconstructor::Visit(int i, int j, int t, Bool3D* not_counted)
    {
        if (i > not_counted->size() - 1 || i < 0)
        {
            return 0;
        }
        else if (j > not_counted->at(i).size() - 1 || j < 0)
        {
            return 0;
        }
        else if (!not_counted->at(i)[j][t])
        {
            return 0;
        }
        else
        {
            not_counted->at(i)[j][t] = false;
            int count = 1;
            count += Visit(i - 1, j - 1, t, not_counted);
            count += Visit(i - 1, j, t, not_counted);
            count += Visit(i - 1, j + 1, t, not_counted);
            count += Visit(i, j - 1, t, not_counted);
            count += Visit(i, j + 1, t, not_counted);
            count += Visit(i + 1, j - 1, t, not_counted);
            count += Visit(i + 1, j, t, not_counted);
            count += Visit(i + 1, j + 1, t, not_counted);
            return count;
        }
    }

    void Reconstructor::BleedTrigger(int i, int j, int t, const Bool3D* three_sigma, Bool3D* good_pixels)
    {
        if (i > three_sigma->size() - 1 || i < 0)
        {
            return;
        }
        else if (j > three_sigma->at(i).size() - 1 || j < 0)
        {
            return;
        }
        else if (t > three_sigma->at(i)[j].size() - 1 || t < 0)
        {
            return;
        }
        else if (!three_sigma->at(i)[j][t])
        {
            return;
        }
        else
        {
            good_pixels->at(i)[j][t] = true;
            BleedTrigger(i - 1, j - 1, t, three_sigma, good_pixels);
            BleedTrigger(i - 1, j, t, three_sigma, good_pixels);
            BleedTrigger(i - 1, j + 1, t, three_sigma, good_pixels);
            BleedTrigger(i, j - 1, t, three_sigma, good_pixels);
            BleedTrigger(i, j + 1, t, three_sigma, good_pixels);
            BleedTrigger(i + 1, j - 1, t, three_sigma, good_pixels);
            BleedTrigger(i + 1, j, t, three_sigma, good_pixels);
            BleedTrigger(i + 1, j + 1, t, three_sigma, good_pixels);
            BleedTrigger(i, j, t - 1, three_sigma, good_pixels);
            BleedTrigger(i, j, t + 1, three_sigma, good_pixels);
        }
    }

    Shower Reconstructor::MakeShower(double t_0, double r_p, double psi, TRotation to_sd_plane)
    {
        // Reconstruct the shower and transform to the world frame (to_sd_plane goes from world frame frame)
        TVector3 shower_direction = to_sd_plane.Inverse() * TVector3(Cos(psi), -Sin(psi), 0);
        if (shower_direction.Z() > 0.0) shower_direction = -shower_direction;
        TVector3 plane_normal = to_sd_plane.Inverse() * TVector3(0, 0, 1);
        TVector3 impact_direction = plane_normal.Cross(shower_direction);
        return Shower(Shower::Params(), r_p * impact_direction, shower_direction, t_0);
    }

    bool Reconstructor::PointWithinView(TVector3 direction, PhotonCount data)
    {
        return direction.Theta() < data.DetectorAxisAngle() - impact_buffer;
    }
}

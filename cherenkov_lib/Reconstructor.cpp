// Reconstructor.cpp
// cherenkov_lib
//
// Author: Matthew Dutson
//
// Implementation of Reconstructor.h

#include <array>
#include <queue>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TMatrixDSymEigen.h>
#include <TFile.h>

#include "Reconstructor.h"

using namespace TMath;

using std::string;
using boost::property_tree::ptree;

namespace cherenkov_simulator
{
    string Reconstructor::Result::Header()
    {
        return "Triggered, " + Shower::Header() + ", Cherenkov, " + Shower::Header();
    }

    string Reconstructor::Result::ToString()
    {
        return Utility::BoolString(trigger) + ", " + mono.ToString() + ", " + Utility::BoolString(impact) +
                ", " + ckv.ToString();
    }

    Reconstructor::Reconstructor(const ptree& config)
    {
        // Construct the ground plane in the world frame.
        TVector3 ground_normal = Utility::ToVector(config.get<string>("surroundings.ground_normal"));
        TVector3 ground_point = Utility::ToVector(config.get<string>("surroundings.ground_point"));
        ground = Plane(ground_normal, ground_point);

        // The rotation from detector frame to world frame. The frames share the same x-axis.
        to_world = Utility::MakeRotation(config.get<double>("surroundings.elevation_angle"));

        // The amount of night sky background noise
        double stop_diameter = config.get<double>("detector.mirror_radius") / (2.0 * config.get<double>("detector.f_number"));
        sky_noise = Sq(stop_diameter / 2.0) * Pi() * global_sky_noise;
        ground_noise = Sq(stop_diameter / 2.0) * Pi() * global_ground_noise;

        // Parameters used when applying triggering logic and noise reduction
        trigger_thresh = config.get<double>("triggering.trigger_thresh");
        noise_thresh = config.get<double>("triggering.noise_thresh");
        trigger_clust = config.get<int>("triggering.trigger_clust");
        impact_buffer = config.get<double>("triggering.impact_buffer");
        plane_dev = config.get<double>("triggering.plane_dev");

        // The random number generator
        rng = TRandom3();
        if (config.get<bool>("simulation.time_seed")) rng.SetSeed();
    }

    Reconstructor::Result Reconstructor::Reconstruct(PhotonCount& data)
    {
        Result result = Result();
        result.trigger = DetectorTriggered(data);
        if (result.trigger)
        {
            TRotation to_sdp = FitSDPlane(data);
            result.mono = MonocularFit(data, to_sdp);
            if (PointWithinView(to_world.Inverse() * result.mono.PlaneImpact(ground), data))
            {
                TVector3 impact;
                if (FindGroundImpact(data, impact))
                {
                    result.ckv = HybridFit(data, impact, to_sdp);
                    result.impact = true;
                }
            }
        }
        return result;
    }

    void Reconstructor::AddNoise(PhotonCount& data)
    {
        PhotonCount::Iterator iter = data.GetIterator();
        while (iter.Next())
        {
            bool toward_ground = ground.InFrontOf(to_world * data.Direction(iter));
            data.AddNoise(toward_ground ? ground_noise : sky_noise, iter, rng);
        }
    }

    void Reconstructor::ClearNoise(PhotonCount& data)
    {
        SubtractAverageNoise(data);
        Bool3D not_visited = GetThresholdMatrices(data, noise_thresh);
        Bool3D triggered = GetThresholdMatrices(data, trigger_thresh);
        Bool3D good_pixels = data.GetFalseMatrix();
        FindPlaneSubset(data, triggered);
        std::vector<bool> trig_state = GetTriggeringState(data);

        std::queue<std::array<ULong, 3>> frontier = std::queue<std::array<ULong, 3>>();
        for (ULong x_trig = 0; x_trig < triggered.size(); x_trig++)
        {
            for (ULong y_trig = 0; y_trig < triggered[x_trig].size(); y_trig++)
            {
                for (ULong t_trig = 0; t_trig < triggered[x_trig][y_trig].size(); t_trig++)
                {
                    if (triggered[x_trig][y_trig][t_trig] && trig_state[t_trig])
                    {
                        frontier.push({x_trig, y_trig, t_trig});
                    }

                    while (!frontier.empty())
                    {
                        std::array<ULong, 3> curr = frontier.front();
                        frontier.pop();
                        ULong x = curr[0];
                        ULong y = curr[1];
                        ULong t = curr[2];
                        not_visited[x][y][t] = false;
                        good_pixels[x][y][t] = true;
                        VisitSpaceAdj(x, y, t, frontier, not_visited);
                        VisitTimeAdj(x, y, t, frontier, not_visited);
                    }
                }
            }
        }
        data.Subset(good_pixels);
    }

    Shower Reconstructor::MonocularFit(const PhotonCount& data, TRotation to_sdp, string graph_file)
    {
        // The functional form of the time profile
        std::stringstream func_string = std::stringstream();
        func_string << "[0] + [1] / (" << Utility::c_cent << ") * tan((pi - [2] - x) / 2)";
        TF1 func = TF1("profile_fit", func_string.str().c_str(), -Pi(), Pi());

        // Set names and initial guesses for parameters
        func.SetParNames("t_0", "r_p", "psi");
        func.SetParameters(0.0, 1e6, PiOver2());

        // Perform the fit and optionally write the fit graph to a file
        TGraphErrors data_graph = GetFitGraph(data, to_sdp);
        data_graph.Fit("profile_fit");
        TF1* result = data_graph.GetFunction("profile_fit");
        if (graph_file != "")
        {
            TFile file(graph_file.c_str(), "RECREATE");
            data_graph.Write("fit_graph");
        }

        // Extract results and make a Shower
        double t_0 = result->GetParameter("t_0");
        double r_p = result->GetParameter("r_p");
        double psi = result->GetParameter("psi");
        return MakeShower(t_0, r_p, psi, to_sdp);
    }

    Shower Reconstructor::HybridFit(const PhotonCount& data, TVector3 impact, TRotation to_sdp, string graph_file)
    {
        // Find the angle of the impact direction with the shower-detector frame x-axis
        double impact_distance = impact.Mag();
        double theta = (to_sdp * impact).Phi();

        // The functional form of the time profile
        std::stringstream func_string = std::stringstream();
        func_string << "[0] +" << impact_distance << " * sin([1] + " << theta << ") / (" << Utility::c_cent << ") * tan((pi - [1] - x) / 2)";
        TF1 func = TF1("profile_fit", func_string.str().c_str(), -Pi(), Pi());

        // Set names and initial guesses for parameters
        func.SetParNames("t_0", "psi");
        func.SetParameters(0.0, PiOver2());

        // Perform the fit and optionally write the fit graph to a file
        TGraphErrors data_graph = GetFitGraph(data, to_sdp);
        data_graph.Fit("profile_fit");
        TF1* result = data_graph.GetFunction("profile_fit");
        if (graph_file != "")
        {
            TFile file(graph_file.c_str(), "RECREATE");
            data_graph.Write("fit_graph");
        }

        // Extract results and make a Shower
        double t_0 = result->GetParameter("t_0");
        double psi = result->GetParameter("psi");
        double r_p = impact_distance * Sin(psi);
        return MakeShower(t_0, r_p, psi, to_sdp);
    }

    TRotation Reconstructor::FitSDPlane(const PhotonCount& data, const Bool3D* mask)
    {
        // Construct the symmetric matrix for finding eigenvectors
        PhotonCount::Iterator iter = data.GetIterator();
        Bool1D curr_mask = Bool1D(data.NBins(), true);
        TMatrixDSym matrix(3);
        for (int j = 0; j < 3; j++)
        {
            for (int k = j; k < 3; k++)
            {
                double mat_element = 0;
                iter.Reset();
                while (iter.Next())
                {
                    int pmt_sum;
                    if (mask == nullptr) pmt_sum = data.SumBins(iter);
                    else pmt_sum = data.FilteredSum(iter, mask->at(iter.X())[iter.Y()]);
                    TVector3 direction = data.Direction(iter);
                    mat_element += direction[j] * direction[k] * pmt_sum;
                }
                matrix[j][k] = mat_element;
                matrix[k][j] = mat_element;
            }
        }

        // Construct a rotation which takes points to the frame where the shower-detector plane is the xy plane.
        TVector3 normal = to_world * MinValVec(matrix);
        if (normal.X() < 0) normal = -normal;
        TVector3 new_x = (normal == TVector3(0, 0, 1)) ? TVector3(1, 0, 0) : TVector3(0, 0, 1).Cross(normal).Unit();
        TVector3 new_y = normal.Cross(new_x).Unit();
        return TRotation().RotateAxes(new_x, new_y, normal).Inverse();
    }

    TVector3 Reconstructor::MinValVec(TMatrixDSym matrix)
    {
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
        return TVector3(eigen_vec[0][min_index], eigen_vec[1][min_index], eigen_vec[2][min_index]);
    }

    bool Reconstructor::FindGroundImpact(PhotonCount& data, TVector3& impact)
    {
        // Find the brightest pixel below the horizon
        TVector3 impact_direction = TVector3();
        int highest_count = 0;
        PhotonCount::Iterator iter = data.GetIterator();
        while (iter.Next())
        {
            TVector3 direction = to_world * data.Direction(iter);
            int sum = data.SumBins(iter);
            if (sum > highest_count && ground.InFrontOf(direction))
            {
                highest_count = sum;
                impact_direction = direction;
            }
        }

        // Find the impact position by extending the direction to the ground
        Ray outward_ray = Ray(TVector3(), impact_direction, 0);
        outward_ray.PropagateToPlane(ground);
        impact = outward_ray.Position();

        // Ensure that the number of photons seen in this brightest pixel exceeds some threshold
        return highest_count > data.FindThreshold(ground_noise, noise_thresh);
    }

    TGraphErrors Reconstructor::GetFitGraph(const PhotonCount& data, TRotation to_sdp)
    {
        // Populate a TGraph with angle/time points from the data. Set errors based on the number of photons viewed.
        Double1D angles = Double1D();
        Double1D times = Double1D();
        Double1D time_er = Double1D();
        PhotonCount::Iterator iter = data.GetIterator();
        while (iter.Next())
        {
            // Don't rotate to the world because the rotation goes from the detector frame to the shower-detector frame.
            TVector3 direction = to_world * data.Direction(iter);
            int bin_sum = data.SumBins(iter);
            if (!ground.InFrontOf(direction) && bin_sum > 0)
            {
                angles.push_back((to_sdp * direction).Phi());
                times.push_back(data.AverageTime(iter));
                time_er.push_back(data.TimeError(iter));
            }
        }

        // Make std::arrays to pass to the TGraph constructor.
        Double1D angle_er = Double1D(angles.size(), 0.0);
        TGraphErrors graph = TGraphErrors((int) angles.size(), &(angles[0]), &(times[0]), &(angle_er[0]), &(time_er[0]));
        graph.Sort();
        return graph;
    }

    void Reconstructor::SubtractAverageNoise(PhotonCount& data)
    {
        PhotonCount::Iterator iter = data.GetIterator();
        while (iter.Next())
        {
            bool toward_ground = ground.InFrontOf(to_world * data.Direction(iter));
            data.Subtract(iter, toward_ground ? ground_noise : sky_noise);
        }
    }

    Bool1D Reconstructor::GetTriggeringState(PhotonCount& data)
    {
        Bool3D trig_matrices = GetThresholdMatrices(data, trigger_thresh, false);
        Bool1D good_frames = Bool1D();
        std::queue<std::array<ULong, 3>> frontier = std::queue<std::array<ULong, 3>>();
        for (ULong x_trig = 0; x_trig < trig_matrices.size(); x_trig++)
        {
            for (ULong y_trig = 0; y_trig < trig_matrices.at(x_trig).size(); y_trig++)
            {
                bool found = false;
                for (ULong t_trig = 0; t_trig < trig_matrices.at(x_trig)[y_trig].size() && !found; t_trig++)
                {
                    int adjacent = 0;
                    if (trig_matrices.at(x_trig)[y_trig][t_trig])
                    {
                        frontier.push({x_trig, y_trig, t_trig});
                    }

                    while (!frontier.empty())
                    {
                        std::array<ULong, 3> curr = frontier.front();
                        ULong x = curr[0];
                        ULong y = curr[1];
                        ULong t = curr[2];
                        trig_matrices.at(x)[y][t] = false;
                        adjacent++;
                        if (adjacent > trigger_clust)
                        {
                            found = true;
                            break;
                        }
                        VisitSpaceAdj(x, y, t, frontier, trig_matrices);
                    }
                }
                good_frames.push_back(found);
            }
        }
        return good_frames;
    }

    void Reconstructor::VisitSpaceAdj(ULong x, ULong y, ULong t, std::queue<std::array<ULong, 3>>& front, const Bool3D& good)
    {
        SafePush(x - 1, y - 1, t, front, good);
        SafePush(x - 1, y, t, front, good);
        SafePush(x - 1, y + 1, t, front, good);
        SafePush(x, y - 1, t, front, good);
        SafePush(x, y + 1, t, front, good);
        SafePush(x + 1, y - 1, t, front, good);
        SafePush(x + 1, y, t, front, good);
        SafePush(x + 1, y + 1, t, front, good);
    }

    void Reconstructor::VisitTimeAdj(ULong x, ULong y, ULong t, std::queue<std::array<ULong, 3>>& front, const Bool3D& good)
    {
        SafePush(x, y, t - 1, front, good);
        SafePush(x, y, t + 1, front, good);
    }

    void Reconstructor::SafePush(ULong x, ULong y, ULong t, std::queue<std::array<ULong, 3>>& front, const Bool3D& good)
    {
        if (x > good.size()) return;
        else if(y > good.at(x).size()) return;
        else if(t > good.at(x)[y].size()) return;
        else if(good.at(x)[y][t]) front.push({x, y, t});
    }

    void Reconstructor::FindPlaneSubset(const PhotonCount& data, Bool3D& triggered)
    {
        TRotation to_sd_plane = FitSDPlane(data, &triggered);
        PhotonCount::Iterator iter = data.GetIterator();
        while (iter.Next())
        {
            if (!NearPlane(to_sd_plane, to_world * data.Direction(iter)))
            {
                triggered.at(iter.X())[iter.Y()] = Bool1D(data.NBins(), false);
            }
        }
    }

    bool Reconstructor::NearPlane(TRotation to_plane, TVector3 direction)
    {
        TVector3 dir_rotated = to_plane * direction;
        double angle = ASin(dir_rotated.Z() / dir_rotated.Mag());
        return Abs(angle) < plane_dev;
    }

    bool Reconstructor::DetectorTriggered(PhotonCount& data)
    {
        Bool1D trig_state = GetTriggeringState(data);
        for (int i = 0; i < trig_state.size(); i++) if (trig_state[i]) return true;
        return false;
    }

    Bool3D Reconstructor::GetThresholdMatrices(PhotonCount& data, double sigma_mult, bool use_below_horiz)
    {
        int ground_thresh = data.FindThreshold(ground_noise, sigma_mult);
        int sky_thresh = data.FindThreshold(sky_noise, sigma_mult);
        Bool3D pass = data.GetFalseMatrix();
        PhotonCount::Iterator iter = data.GetIterator();
        while (iter.Next())
        {
            bool toward_ground = ground.InFrontOf(to_world * data.Direction(iter));
            if (toward_ground && !use_below_horiz) continue;
            pass[iter.X()][iter.Y()] = data.AboveThreshold(iter, toward_ground ? ground_thresh : sky_thresh);
        }
        return pass;
    }

    Shower Reconstructor::MakeShower(double t_0, double r_p, double psi, TRotation to_sd_plane)
    {
        // Reconstruct the shower and transform to the world frame (to_sd_plane goes from world frame)
        TVector3 shower_direction = to_sd_plane.Inverse() * TVector3(Cos(psi), -Sin(psi), 0);
        if (shower_direction.Z() > 0.0) shower_direction = -shower_direction;
        TVector3 plane_normal = to_sd_plane.Inverse() * TVector3(0, 0, 1);
        TVector3 impact_direction = plane_normal.Cross(shower_direction);
        return Shower(Shower::Params(), r_p * impact_direction, shower_direction, t_0);
    }

    bool Reconstructor::PointWithinView(TVector3 direction, const PhotonCount& data)
    {
        return direction.Theta() < data.DetectorAxisAngle() - impact_buffer;
    }
}

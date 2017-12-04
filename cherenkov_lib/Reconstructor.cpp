// Reconstructor.cpp
//
// Author: Matthew Dutson
//
// Implementation of Reconstructor.h

#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TMatrixDSymEigen.h>

#include "Reconstructor.h"

using namespace std;
using namespace boost::property_tree;
using namespace TMath;

namespace cherenkov_simulator
{
    Reconstructor::Result::Result()
    {
        triggered = false;
        chkv_tried = false;
    }

    string Reconstructor::Result::Header()
    {
        return "Triggered," + Shower::Header() + ",Cherenkov," + Shower::Header();
    }

    string Reconstructor::Result::ToString(Plane ground_plane) const
    {
        string result;
        if (triggered) result += "1," + mono_recon.ToString(ground_plane) + ",";
        else result += "0,0,0,0,";

        if (chkv_tried) result += "1," + chkv_recon.ToString(ground_plane);
        else result += "0,0,0,0";
        return result;
    }

    Reconstructor::Reconstructor(const ptree& config)
    {
        TVector3 ground_norm = Utility::ToVector(config.get<string>("surroundings.ground_norm"));
        TVector3 ground_fixd = Utility::ToVector(config.get<string>("surroundings.ground_fixd"));
        ground_plane = Plane(ground_norm, ground_fixd);
        rot_to_world = Utility::MakeRotation(config.get<double>("surroundings.elevation_angle"));

        double mirror_radius = config.get<double>("detector.mirror_radius");
        double stop_diameter = mirror_radius / (2.0 * config.get<double>("detector.f_number"));
        sky_noise = Sq(stop_diameter / 2.0) * Pi() * glob_sky_noise;
        gnd_noise = Sq(stop_diameter / 2.0) * Pi() * glob_gnd_noise;

        trigr_thresh = config.get<double>("triggering.trigr_thresh");
        noise_thresh = config.get<double>("triggering.noise_thresh");
        impact_buffr = config.get<double>("triggering.impact_buffr");
        plane_thresh = config.get<double>("triggering.plane_thresh");
        trigr_clustr = config.get<int>("triggering.trigr_clustr");
    }

    Reconstructor::Result Reconstructor::Reconstruct(const PhotonCount& data) const
    {
        Result result = Result();
        result.triggered = DetectorTriggered(GetTriggeringState(data));
        if (result.triggered)
        {
            TRotation to_sdp = FitSDPlane(data);
            result.mono_recon = MonocularFit(data, to_sdp);
            TVector3 direction = rot_to_world.Inverse() * result.mono_recon.PlaneImpact(ground_plane);
            if (direction.Theta() < data.DetectorAxisAngle() - impact_buffr)
            {
                TVector3 impact;
                if (FindGroundImpact(data, impact))
                {
                    result.chkv_recon = HybridFit(data, impact, to_sdp);
                    result.chkv_tried = true;
                }
            }
        }
        return result;
    }

    void Reconstructor::AddNoise(PhotonCount& data) const
    {
        PhotonCount::Iterator iter = data.GetIterator();
        while (iter.Next())
        {
            bool toward_ground = ground_plane.InFrontOf(rot_to_world * data.Direction(iter));
            data.AddNoise(toward_ground ? gnd_noise : sky_noise, iter);
        }
    }

    Shower Reconstructor::MonocularFit(const PhotonCount& data, TRotation to_sdp, string graph_file) const
    {
        // Define the functional form of the time profile.
        stringstream func_string = stringstream();
        func_string << "[0] + [1] / (" << c_cent << ") * tan((pi - [2] - x) / 2)";
        TF1 func = TF1("profile_fit", func_string.str().c_str(), -Pi(), Pi());

        func.SetParNames("t_0", "r_p", "psi");
        func.SetParameters(0.0, 1e6, PiOver2());
        TGraphErrors data_graph = GetFitGraph(data, to_sdp);
        data_graph.Fit("profile_fit", "Q");
        TF1* result = data_graph.GetFunction("profile_fit");
        if (!graph_file.empty())
        {
            TFile file(graph_file.c_str(), "RECREATE");
            data_graph.Write("fit_graph");
        }

        double t_0 = result->GetParameter("t_0");
        double r_p = result->GetParameter("r_p");
        double psi = result->GetParameter("psi");
        return MakeShower(t_0, r_p, psi, to_sdp);
    }

    Shower Reconstructor::HybridFit(const PhotonCount& data, TVector3 impact, TRotation to_sdp, string graph_file) const
    {
        double impact_distance = impact.Mag();
        double alpha = (to_sdp * impact).Phi();

        // Define the functional form of the time profile.
        stringstream func_string = stringstream();
        func_string << "[0] +" << impact_distance << " * sin([1] + " << alpha << ") / (" << c_cent
                    << ") * tan((pi - [1] - x) / 2)";
        TF1 func = TF1("profile_fit", func_string.str().c_str(), -Pi(), Pi());

        func.SetParNames("t_0", "psi");
        func.SetParameters(0.0, PiOver2());
        TGraphErrors data_graph = GetFitGraph(data, to_sdp);
        data_graph.Fit("profile_fit", "Q");
        TF1* result = data_graph.GetFunction("profile_fit");
        if (!graph_file.empty())
        {
            TFile file(graph_file.c_str(), "RECREATE");
            data_graph.Write("fit_graph");
        }

        double t_0 = result->GetParameter("t_0");
        double psi = result->GetParameter("psi");
        double r_p = impact_distance * Sin(psi);
        return MakeShower(t_0, r_p, psi, to_sdp);
    }

    TRotation Reconstructor::FitSDPlane(const PhotonCount& data, const Bool3D* mask) const
    {
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
                    if (mask == nullptr)
                        pmt_sum = data.SumBins(iter);
                    else
                        pmt_sum = data.SumBinsFiltered(iter, *mask);
                    TVector3 direction = data.Direction(iter);
                    mat_element += direction[j] * direction[k] * pmt_sum;
                }
                matrix[j][k] = mat_element;
                matrix[k][j] = mat_element;
            }
        }

        // Construct the shower-detector frame.
        TVector3 normal = rot_to_world * MinValVec(matrix);
        if (normal.X() < 0) normal = -normal;
        TVector3 new_x = (normal == TVector3(0, 0, 1)) ? TVector3(1, 0, 0) : TVector3(0, 0, 1).Cross(normal).Unit();
        TVector3 new_y = normal.Cross(new_x).Unit();
        return TRotation().RotateAxes(new_x, new_y, normal).Inverse();
    }

    TVector3 Reconstructor::MinValVec(TMatrixDSym matrix) const
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

    bool Reconstructor::FindGroundImpact(const PhotonCount& data, TVector3& impact) const
    {
        TVector3 reflect_dir = TVector3();
        int highest_sum = 0;
        PhotonCount::Iterator iter = data.GetIterator();
        while (iter.Next())
        {
            TVector3 direction = rot_to_world * data.Direction(iter);
            int sum = data.SumBins(iter);
            if (sum > highest_sum && ground_plane.InFrontOf(direction))
            {
                highest_sum = sum;
                reflect_dir = direction;
            }
        }

        if (reflect_dir == TVector3()) return false;
        Ray outward_ray = Ray(TVector3(), reflect_dir, 0);
        outward_ray.PropagateToPlane(ground_plane);
        impact = outward_ray.Position();
        return highest_sum > data.FindThreshold(gnd_noise, trigr_thresh);
    }

    TGraphErrors Reconstructor::GetFitGraph(const PhotonCount& data, TRotation to_sdp) const
    {
        Double1D angles = Double1D();
        Double1D times = Double1D();
        Double1D time_err = Double1D();
        PhotonCount::Iterator iter = data.GetIterator();
        while (iter.Next())
        {
            // Don't rotate to the world because the rotation goes from the detector frame to the shower-detector frame.
            TVector3 direction = rot_to_world * data.Direction(iter);
            int bin_sum = data.SumBins(iter);
            if (!ground_plane.InFrontOf(direction) && bin_sum > 0)
            {
                angles.push_back((to_sdp * direction).Phi());
                times.push_back(data.AverageTime(iter));
                time_err.push_back(data.TimeError(iter));
            }
        }

        Double1D angle_err = Double1D(angles.size(), 0.0);
        TGraphErrors graph = TGraphErrors((int) angles.size(), &(angles[0]), &(times[0]), &(angle_err[0]), &(time_err[0]));
        graph.Sort();
        return graph;
    }

    void Reconstructor::SubtractAverageNoise(PhotonCount& data) const
    {
        PhotonCount::Iterator iter = data.GetIterator();
        while (iter.Next())
        {
            bool toward_ground = ground_plane.InFrontOf(rot_to_world * data.Direction(iter));
            data.Subtract(toward_ground ? gnd_noise : sky_noise, iter);
        }
    }

    Bool1D Reconstructor::GetTriggeringState(const PhotonCount& data) const
    {
        Bool3D trig_matrices = GetThresholdMatrices(data, trigr_thresh, false);
        Bool1D good_frames = Bool1D(data.NBins(), false);

        list<array<size_t, 3>> frontier = list<array<size_t, 3>>();
        for (size_t t_trig = 0; t_trig < data.NBins(); t_trig++)
        {
            bool found = false;
            for (size_t x_trig = 0; x_trig < trig_matrices.size() && !found; x_trig++)
            {
                for (size_t y_trig = 0; y_trig < trig_matrices.at(x_trig).size() && !found; y_trig++)
                {
                    if (trig_matrices.at(x_trig)[y_trig][t_trig])
                        frontier.push_back({x_trig, y_trig, t_trig});

                    int adjacent = 0;
                    while (!frontier.empty() && !found)
                    {
                        array<size_t, 3> curr = frontier.front();
                        frontier.pop_front();
                        size_t x = curr[0];
                        size_t y = curr[1];
                        size_t t = curr[2];
                        adjacent++;
                        if (adjacent > trigr_clustr)
                        {
                            frontier.clear();
                            found = true;
                        }
                        VisitSpaceAdj(x, y, t, frontier, trig_matrices);
                    }
                }
            }
            good_frames[t_trig] = found;
        }
        return good_frames;
    }

    void Reconstructor::ClearNoise(PhotonCount& data) const
    {
        SubtractAverageNoise(data);
        Bool3D not_visited = GetThresholdMatrices(data, noise_thresh);
        Bool3D triggered = GetThresholdMatrices(data, trigr_thresh);
        Bool3D good_pixels = data.GetFalseMatrix();
        FindPlaneSubset(data, triggered);
        Bool1D trig_state = GetTriggeringState(data);

        list<array<size_t, 3>> frontier = list<array<size_t, 3>>();
        for (size_t x_trig = 0; x_trig < triggered.size(); x_trig++)
        {
            for (size_t y_trig = 0; y_trig < triggered[x_trig].size(); y_trig++)
            {
                for (size_t t_trig = 0; t_trig < triggered[x_trig][y_trig].size(); t_trig++)
                {
                    if (triggered[x_trig][y_trig][t_trig] && trig_state[t_trig])
                        frontier.push_back({x_trig, y_trig, t_trig});

                    while (!frontier.empty())
                    {
                        array<size_t, 3> curr = frontier.front();
                        frontier.pop_front();
                        size_t x = curr[0];
                        size_t y = curr[1];
                        size_t t = curr[2];
                        good_pixels[x][y][t] = true;
                        VisitSpaceAdj(x, y, t, frontier, not_visited);
                        VisitTimeAdj(x, y, t, frontier, not_visited);
                    }
                }
            }
        }
        data.Subset(good_pixels);
    }

    void Reconstructor::VisitSpaceAdj(size_t x, size_t y, size_t t, list<array<size_t, 3>>& front, Bool3D& not_visited)
    {
        VisitPush(x - 1, y - 1, t, front, not_visited);
        VisitPush(x - 1, y, t, front, not_visited);
        VisitPush(x - 1, y + 1, t, front, not_visited);
        VisitPush(x, y - 1, t, front, not_visited);
        VisitPush(x, y + 1, t, front, not_visited);
        VisitPush(x + 1, y - 1, t, front, not_visited);
        VisitPush(x + 1, y, t, front, not_visited);
        VisitPush(x + 1, y + 1, t, front, not_visited);
    }

    void Reconstructor::VisitTimeAdj(size_t x, size_t y, size_t t, list<array<size_t, 3>>& front, Bool3D& not_visited)
    {
        VisitPush(x, y, t - 1, front, not_visited);
        VisitPush(x, y, t + 1, front, not_visited);
    }

    void Reconstructor::VisitPush(size_t x, size_t y, size_t t, list<array<size_t, 3>>& front, Bool3D& not_visited)
    {
        if (x > not_visited.size() - 1)
            return;
        else if (y > not_visited.at(x).size() - 1)
            return;
        else if (t > not_visited.at(x)[y].size() - 1)
            return;
        else if (not_visited.at(x)[y][t])
            front.push_back({x, y, t});
        not_visited[x][y][t] = false;
    }

    void Reconstructor::FindPlaneSubset(const PhotonCount& data, Bool3D& triggered) const
    {
        TRotation to_sd_plane = FitSDPlane(data, &triggered);
        PhotonCount::Iterator iter = data.GetIterator();
        while (iter.Next())
        {
            if (!NearPlane(to_sd_plane, rot_to_world * data.Direction(iter)))
                triggered[iter.X()][iter.Y()] = Bool1D(data.NBins(), false);
        }
    }

    bool Reconstructor::NearPlane(TRotation to_plane, TVector3 direction) const
    {
        TVector3 dir_rotated = to_plane * direction;
        double angle = ASin(dir_rotated.Z() / dir_rotated.Mag());
        return Abs(angle) < plane_thresh;
    }

    bool Reconstructor::DetectorTriggered(const Bool1D& trig_state) const
    {
        for (bool state : trig_state)
            if (state) return true;
        return false;
    }

    Bool3D Reconstructor::GetThresholdMatrices(const PhotonCount& data, double sigma_mult, bool use_below_horiz) const
    {
        int gnd_thresh = data.FindThreshold(gnd_noise, sigma_mult);
        int sky_thresh = data.FindThreshold(sky_noise, sigma_mult);
        Bool3D pass = data.GetFalseMatrix();
        PhotonCount::Iterator iter = data.GetIterator();
        while (iter.Next())
        {
            bool toward_ground = ground_plane.InFrontOf(rot_to_world * data.Direction(iter));
            if (toward_ground && !use_below_horiz) continue;
            pass[iter.X()][iter.Y()] = data.AboveThreshold(iter, toward_ground ? gnd_thresh : sky_thresh);
        }
        return pass;
    }

    Shower Reconstructor::MakeShower(double t_0, double r_p, double psi, TRotation to_sdp)
    {
        // Remember that to_sdp goes from the world frame.
        TVector3 shower_direction = to_sdp.Inverse() * TVector3(Cos(psi), -Sin(psi), 0);
        if (shower_direction.Z() > 0.0) shower_direction = -shower_direction;
        TVector3 plane_normal = to_sdp.Inverse() * TVector3(0, 0, 1);
        TVector3 impact_direction = plane_normal.Cross(shower_direction);
        return Shower(1.0, 1.0, r_p * impact_direction, shower_direction, t_0);
    }
}

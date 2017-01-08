// common.cpp
// cherenkov_simulator
//
// Created by Matthew Dutson on 9/8/16.
//
//

#include "utility.h"

using namespace boost::program_options;
using namespace std;

namespace cherenkov_simulator
{
    TVector3 ToVector(string s)
    {
        TVector3 output = TVector3();
        int comma = s.find(',');
        output.SetX(stod(s.substr(1, comma)));
        s.erase(0, comma);
        comma = s.find(',');
        output.SetY(stod(s.substr(0, comma)));
        s.erase(0, comma);
        output.SetZ(stod(s.substr(0, s.size() - 1)));
        return output;
    }

    TVector3 RandomPerpendicularVector(TVector3 vec, TRandom3 rng)
    {
        if (vec.X() == 0 && vec.Y() == 0 && vec.Z() == 0)
        {
            return TVector3(1, 0, 0);
        } else
        {
            TVector3 other_vec = vec + TVector3(1, 0, 0);
            TVector3 normal = (vec.Cross(other_vec)).Unit();
            normal.Rotate(rng.Uniform(2 * TMath::Pi()), vec);
            return normal;
        }
    }
}

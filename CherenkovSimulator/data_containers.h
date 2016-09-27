// data_containers.h
// cherenkov_simulator
//
// Created by Matthew Dutson on 9/8/16.
//
//

#ifndef data_containers_h
#define data_containers_h

namespace cherenkov_simulator
{
    class RawData
    {
    private:
        
        std::vector<TVector3> directions;
        
        std::vector<double> times;
        
    public:
        
        void AddPoint(TVector3 direction, double time);
    };

    class SignalData
    {
        
    };
}

#endif

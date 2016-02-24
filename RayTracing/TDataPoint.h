/*
 * Created by Matthew Dutson on 2/23/16.
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * This class represents a data point with a time, an x position, and a y position.
 */

#ifndef TDataPoint_h
#define TDataPoint_h

#include "TMath.h"

class TDataPoint {
    
private:
    
    // The time at which the data point was collected
    Double_t fTime;
    
    // The x position
    Double_t fX;
    
    // The y position
    Double_t fY;
    
public:
    
    /*
     * The standard constructor.
     */
    TDataPoint(Double_t time, Double_t x, Double_t y);
    
    /*
     * Returns the time stored in the data point.
     */
    Double_t GetTime();
    
    /*
     * Returns the x position of the data point.
     */
    Double_t GetX();
    
    /*
     * Returns the y position of the data point.
     */
    Double_t GetY();
    
};

#endif

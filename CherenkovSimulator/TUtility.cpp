//
//  TUtility.cpp
//  CherenkovSimulator
//
//  Created by Matthew Dutson on 3/22/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#include "TUtility.h"

Double_t TUtility::SumArray(std::vector<Double_t> array) {
    Double_t sum = 0;
    for (Double_t d: array) {
        sum += d;
    }
    return sum;
}
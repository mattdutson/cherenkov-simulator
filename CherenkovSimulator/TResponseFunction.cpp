/*
 * CherenkovSimulator - TCoordinates.cpp
 *
 * Copyright © 2016 Matthew Dutson. All rights reserved.
 *
 * Contains the implementation of TResponseFunction.
 */

#include "TResponseFunction.hpp"

TResponseFunction::TResponseFunction(): TF1("default", "1") {
    fResponseTime = 1;
}

TResponseFunction::TResponseFunction(Double_t responseTime, TF1 responseFunction): TF1(responseFunction) {
    fResponseTime = responseTime;
}

TH1D TResponseFunction::ResponseHistogram(Double_t sampleTime) {
    TH1D output = TH1D("Response Histogram", "Response Histogram", fResponseTime / sampleTime, 0, fResponseTime);
    for (Int_t i = 0; i < output.GetSize(); i++) {
        Double_t t = i * sampleTime;
        output.SetBinContent(i, Eval(t));
    }
    return output;
}

Double_t TResponseFunction::ResponseTime() {
    return fResponseTime;
}
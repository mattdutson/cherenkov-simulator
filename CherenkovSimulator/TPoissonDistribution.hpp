//
//  TPoissonDistribution.hpp
//  CherenkovSimulator
//
//  Created by Matthew Dutson on 4/13/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#ifndef TPoissonDistribution_hpp
#define TPoissonDistribution_hpp

#include <stdio.h>
#include "TF1.h"
#include "TRandom1.h"

class TPoissonDistribution: TF1 {
    
private:
    
    Double_t fAverageNumber;
    
    std::vector<Double_t> fIntegratedDistro;
    
    TRandom1* fRng = new TRandom1(12342834, 3);
    
public:
    
    TPoissonDistribution(Double_t averageNumber);
    
    Int_t GetNumber();
};

#endif /* TPoissonDistribution_hpp */

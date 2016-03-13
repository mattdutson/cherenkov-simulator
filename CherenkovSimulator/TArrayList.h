//
//  TArrayList.hpp
//  RayTracing
//
//  Created by Matthew Dutson on 2/10/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#ifndef TArrayList_h
#define TArrayList_h

#include <stdio.h>
#include "TMath.h"

template <class Type> class TArrayList {
private:
    
    Int_t size;
    Type contents[];
    
    void DoubleArraySize();
    
public:
    
    TArrayList();
    
    void AddLast(Type item);
    
    Int_t GetSize();
    
    Type * GetContents();
};

#endif /* TArrayList_hpp */

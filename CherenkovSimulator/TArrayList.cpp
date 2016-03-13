//
//  TArrayList.cpp
//  RayTracing
//
//  Created by Matthew Dutson on 2/10/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#include "TArrayList.h"

template <class Type> TArrayList<Type>::TArrayList() {
    contents = Type[128];
    
}
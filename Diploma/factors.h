//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#ifndef __Diploma__factors__
#define __Diploma__factors__

#include <iostream>

namespace ftr {

    double cEf(double T, double dT);
    
    double alpha(double t);
    double sigma(double t);

    extern double const X1;
    extern double const X2;
    extern double const totalTime;

    extern double const TStart;
    extern double const TEnv;
    extern double const TEnv4;

};

#endif /* defined(__Diploma__factors__) */

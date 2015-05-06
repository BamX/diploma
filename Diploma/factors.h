//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#ifndef __Diploma__factors__
#define __Diploma__factors__

#include "config.h"

class Factors {
    Config config;
public:

    double cEf(double T) const;
    
    double alpha(double t) const;
    double sigma(double t) const;

    double lambda(double T) const;
    double ro(double T) const;

    double X1() const;
    double X2() const;
    double totalTime() const;

    double TStart() const;
    double TEnv() const;
    double TEnv4() const;

    size_t ViewCount() const;
    size_t DebugView() const;
    size_t MatrixFramesCount() const;
    double X1View(size_t index) const;
    double X2View(size_t index) const;

};

#endif /* defined(__Diploma__factors__) */

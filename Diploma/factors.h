//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#ifndef __Diploma__factors__
#define __Diploma__factors__

#include "config.h"
#include <vector>

class Factors {
    Config config;

    double _x1, _x2, _totalTime,
        _x1SplitCount, _x2SplitCount, _timeSplitCount, _epsilon, _tMax,
        _TStart, _TEnv, _TEnv4;
    bool _enableConsole, _enablePlot, _enableMatrix;
    size_t _viewCount, _debugView, _framesCount;
    std::vector<double> _x1View, _x2View;

public:

    Factors();

    double cEf(double T) const;
    
    double alpha(double t) const;
    double sigma(double t) const;

    double lambda(double T) const;
    double ro(double T) const;

    double X1() const;
    double X2() const;
    double totalTime() const;

    double X1SplitCount() const;
    double X2SplitCount() const;
    double TimeSplitCount() const;
    double Epsilon() const;
    double TMax() const;

    double TStart() const;
    double TEnv() const;
    double TEnv4() const;

    bool EnableConsole() const;
    bool EnablePlot() const;
    bool EnableMatrix() const;

    size_t ViewCount() const;
    size_t DebugView() const;
    size_t FramesCount() const;
    double X1View(size_t index) const;
    double X2View(size_t index) const;

};

#endif /* defined(__Diploma__factors__) */

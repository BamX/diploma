//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#ifndef __Diploma__factors__
#define __Diploma__factors__

#include "config.h"
#include <vector>

class Factors {
    Config config;

    float _x1, _x2, _totalTime,
        _x1SplitCount, _x2SplitCount, _timeSplitCount, _epsilon,
        _TStart, _TEnv, _TEnv4;
    bool _enableConsole, _enablePlot, _enableMatrix;
    size_t _viewCount, _debugView, _matrixFramesCount;
    std::vector<float> _x1View, _x2View;

public:

    Factors();

    float cEf(float T) const;
    
    float alpha(float t) const;
    float sigma(float t) const;

    float lambda(float T) const;
    float ro(float T) const;

    float X1() const;
    float X2() const;
    float totalTime() const;

    float X1SplitCount() const;
    float X2SplitCount() const;
    float TimeSplitCount() const;
    float Epsilon() const;

    float TStart() const;
    float TEnv() const;
    float TEnv4() const;

    bool EnableConsole() const;
    bool EnablePlot() const;
    bool EnableMatrix() const;

    size_t ViewCount() const;
    size_t DebugView() const;
    size_t MatrixFramesCount() const;
    float X1View(size_t index) const;
    float X2View(size_t index) const;

};

#endif /* defined(__Diploma__factors__) */

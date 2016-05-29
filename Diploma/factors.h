//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#ifndef __Diploma__factors__
#define __Diploma__factors__

#include "config.h"
#include <vector>

class Factors {
    Config *_config;

    double _x1, _x2, _totalTime,
        _x1SplitCount, _x2SplitCount, _timeSplitCount, _epsilon, _tMax,
        _TStart, _TEnv, _TEnv4, _balanceFactor;
    bool _enableConsole, _enablePlot, _enableMatrix, _enableBuckets;
    size_t _minimumBundle, _viewCount, _debugView, _framesCount, _repeats;
    std::vector<double> _x1View, _x2View;
    std::string _plotFilename, _bucketsFilename;

    void initFactors(Config config);

public:

    ~Factors();
    void initFactors(const char *filename);

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
    size_t Repeats() const;

    size_t MinimumBundle() const;
    double BalanceFactor() const;

    double TStart() const;
    double TEnv() const;
    double TEnv4() const;

    bool EnableConsole() const;
    bool EnablePlot() const;
    bool EnableMatrix() const;
    bool EnableBuckets() const;

    std::string PlotFilename() const;
    std::string MatrixFilename() const;
    std::string BucketsFilename() const;

    size_t ViewCount() const;
    size_t DebugView() const;
    size_t FramesCount() const;
    double X1View(size_t index) const;
    double X2View(size_t index) const;

};

#endif /* defined(__Diploma__factors__) */

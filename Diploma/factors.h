//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#ifndef __Diploma__factors__
#define __Diploma__factors__

#include "config.h"
#include <vector>

extern size_t const kAlgorithmTranspose;
extern size_t const kAlgorithmStatic;

class Factors {
    Config *_config;

    double _x1, _x2, _totalTime,
        _x1SplitCount, _x2SplitCount, _timeSplitCount, _epsilon, _tMax,
        _TStart, _TEnv, _TEnv4, _balanceFactor, _transposeBalancingFactor,
        _transposeBalancingTimeFactor, _staticBalancingThresholdFactor;
    bool _balancing, _enableConsole, _enablePlot, _enableMatrix, _enableBuckets, _enableWeights, _enableTimes,
        _enableBalanceWeightsSmooth;
    size_t _minimumBundle, _viewCount, _debugView, _framesCount, _repeats, _transposeIterations, _algorithm;
    std::vector<double> _x1View, _x2View;
    std::string _plotFilename, _bucketsFilename, _weightsFilename, _timesFilenamePrefix;

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
    bool Balancing() const;
    double TransposeBalanceFactor() const;
    size_t TransposeBalanceIterationsInterval() const;
    double TransposeBalanceTimeFactor() const;
    double StaticBalanceThresholdFactor() const;
    bool EnableBalanceWeightsSmooth() const;

    size_t Algorithm() const;

    double TStart() const;
    double TEnv() const;
    double TEnv4() const;

    bool EnableConsole() const;
    bool EnablePlot() const;
    bool EnableMatrix() const;
    bool EnableBuckets() const;
    bool EnableWeights() const;
    bool EnableTimes() const;

    std::string PlotFilename() const;
    std::string MatrixFilename() const;
    std::string BucketsFilename() const;
    std::string WeightsFilename() const;
    std::string TimesFilenamePrefix() const;

    size_t ViewCount() const;
    size_t DebugView() const;
    size_t FramesCount() const;
    double X1View(size_t index) const;
    double X2View(size_t index) const;

};

#endif /* defined(__Diploma__factors__) */

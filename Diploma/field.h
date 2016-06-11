//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#ifndef __Diploma__field__
#define __Diploma__field__

#include <fstream>
#include <mpi.h>
#include <vector>
#include <chrono>

extern int const MASTER;
extern int const WAITER;
extern int const NOBODY;
extern int const NOTHING;
extern size_t const MAX_ITTERATIONS_COUNT;

typedef std::chrono::high_resolution_clock bx_clock_t;
typedef bx_clock_t::time_point bx_time_t;

class Field {
protected:
    std::ofstream *fout, *mfout, *bfout, *wfout, *tfout;
    bx_time_t startSyncTime;
    double fullCalculationTime;
    double nextFrameTime;

    double *prev, *curr, *buff, *views;
    double *maF, *mbF, *mcF, *mfF;

    size_t width, height, origWidth, origHeight;
    bool transposed;

    double t;
    double hX, hY, dT;
    double epsilon;
    size_t lastIterrationsCount;

    int myId, numProcs, myCoord;
    MPI_Comm comm;
    size_t mySX, mySY;
    int topN, bottomN, leftN, rightN;

    void fillInitial();
    virtual void calculateNBS();

    void fillFactors(size_t row, bool first);
    void firstPass(size_t row);
    double secondPass(size_t row, bool first);
    double solve(size_t row, bool first);

    virtual size_t solveRows();

    virtual void transpose();
    void nextTimeLayer();

    void enablePlotOutput();
    void enableMatrixOutput();
    void enableBucketsOutput();
    void enableWeightsOutput();
    void enableTimesOutput();

    void printAll();
    virtual void printConsole();
    virtual void printMatrix();
    void printViews();

    unsigned long long picosecFromStart();
    std::ostream &debug(bool info = true);

    void reduceViews();

#pragma mark - Balancing MPI

    double *weights;

    virtual void syncWeights();
    virtual bool balanceNeeded();
    virtual void balance();

    void smoothWeights();

    virtual bool isBucketsMaster();
    virtual size_t weightsSize();

public:
    Field();
    ~Field();
    virtual void finalize();

    void test();
    void testPrint();

    virtual void init();
    void solve();
    double time();
    bool done();

    double calculationTime();

    virtual double view(double x1, double x2);
    double view(size_t index);

};

#endif /* defined(__Diploma__field__) */

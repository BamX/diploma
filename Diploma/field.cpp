//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#include "field.h"
#include "algo.h"
#include <cmath>

int const MASTER = 0;
int const WAITER = 0;
int const NOBODY = MPI_PROC_NULL;
int const NOTHING = -1;

size_t const MAX_ITTERATIONS_COUNT = 50;

Field::Field() {
    fout = NULL;
    mfout = NULL;
    bfout = NULL;
}

Field::~Field() {
    if (fout != NULL) {
        fout->close();
        delete fout;
    }

    if (mfout != NULL) {
        mfout->close();
        delete mfout;
    }
    
    delete[] prev;
    delete[] curr;
    delete[] buff;
    delete[] views;

    delete[] maF;
    delete[] mbF;
    delete[] mcF;
    delete[] mfF;
}

void Field::init() {
    width = algo::ftr().X1SplitCount();
    height = algo::ftr().X2SplitCount();

    hX = algo::ftr().X1() / (width - 1);
    hY = algo::ftr().X2() / (height - 1);
    dT = algo::ftr().totalTime() / algo::ftr().TimeSplitCount();
    epsilon = algo::ftr().Epsilon();
    transposed = false;
    fout = NULL;
    mfout = NULL;
    bfout = NULL;
    wfout = NULL;
    tfout = NULL;

    calculateNBS();

    prev = new double[width * width];
    curr = new double[width * width];
    buff = new double[width * width];
    views = new double[algo::ftr().ViewCount()];

    maF = new double[width * width];
    mbF = new double[width * width];
    mcF = new double[width * width];
    mfF = new double[width * width];

    fillInitial();

    if (algo::ftr().EnablePlot()) {
        enablePlotOutput();
    }
    if (algo::ftr().EnableMatrix()) {
        enableMatrixOutput();
    }
    if (algo::ftr().EnableBuckets()) {
        enableBucketsOutput();
    }
    if (algo::ftr().EnableWeights()) {
        enableWeightsOutput();
    }
    if (algo::ftr().EnableTimes()) {
        enableTimesOutput();
    }

    debug(0).flush();
}

void Field::fillInitial() {
    t = 0;
    nextFrameTime = 0;
    lastIterrationsCount = 0;
    fullProcessingTime = 0;
    if (transposed) {
        transpose();
    }

    for (size_t index = 0, len = width * height; index < len; ++index) {
        curr[index] = algo::ftr().TStart();
    }
}

void Field::transpose() {
}

void Field::nextTimeLayer() {
    std::swap(curr, prev);
    t += dT;
}

void Field::fillFactors(size_t row, bool first) {
    START_TIME(start);
    
    double *aF = maF + row * width;
    double *bF = mbF + row * width;
    double *cF = mcF + row * width;
    double *fF = mfF + row * width;

    double *rw = prev + row * width;
    double *brw = first ? rw : (curr + row * width);

    algo::fillFactors(rw, brw, width, aF, bF, cF, fF, t, hX, dT, leftN == NOBODY, rightN == NOBODY);

    END_TIME(calculationsTime, start);
}

void Field::firstPass(size_t row) {
    START_TIME(start);

    double *aF = maF + row * width;
    double *bF = mbF + row * width;
    double *cF = mcF + row * width;
    double *fF = mfF + row * width;

    algo::firstPass(width, aF, bF, cF, fF, rightN == NOBODY);

    END_TIME(calculationsTime, start);
}

double Field::secondPass(size_t row, bool first) {
    START_TIME(start);

    double *bF = mbF + row * width;
    double *cF = mcF + row * width;
    double *fF = mfF + row * width;

    double *y = curr + row * width;
    double *py = first ? (prev + row * width) : y;

    double maxDelta = 0;
    algo::secondPass(py, y, width, bF, cF, fF, rightN == NOBODY, &maxDelta);

    END_TIME(calculationsTime, start);

    return maxDelta;
}

double Field::solve(size_t row, bool first) {
    firstPass(row);
    return secondPass(row, first);
}

size_t Field::solveRows() {
    return 0;
}

void Field::solve() {
    if (done()) {
        return;
    }

    START_TIME(solveStart);

    lastIterrationsCount = 0;

    nextTimeLayer();
    lastIterrationsCount += solveRows();
    if (balanceNeeded()) {
        syncWeights();
        balance();
    }

    nextTimeLayer();
    transpose();
    lastIterrationsCount += solveRows();
    if (balanceNeeded()) {
        syncWeights();
        balance();
    }

    transpose();

    END_TIME(fullIterationTime, solveStart);

    printAll();

    END_TIME(fullProcessingTime, solveStart) * 1e-12;
}

double Field::time() {
    return t;
}

bool Field::done() {
    return t >= algo::ftr().TMax();
}

double Field::calculationTime()
{
    return fullProcessingTime;
}

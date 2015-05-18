//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#include "field.h"
#include "factors.h"
#include <cmath>

int const MASTER = 0;
int const WAITER = 3;
int const NOBODY = MPI_PROC_NULL;
int const NOTHING = -1;
int const SEND_PACK_SIZE = 6;

static size_t const MAX_ITTERATIONS_COUNT = 50;

Field::Field() {
    initFactors();

    width = ftr.X1SplitCount();
    height = ftr.X2SplitCount();

    hX = ftr.X1() / (width - 1);
    hY = ftr.X2() / (height - 1);
    dT = ftr.totalTime() / ftr.TimeSplitCount();
    epsilon = ftr.Epsilon();
    transposed = false;
    fout = NULL;
    mfout = NULL;

    calculateNBS();

    prev = new double[height * width];
    curr = new double[height * width];
    buff = new double[height * width];
    views = new double[ftr.ViewCount()];

    maF = new double[height * width];
    mbF = new double[height * width];
    mcF = new double[height * width];
    mfF = new double[height * width];

    if (ftr.EnablePlot()) {
        enablePlotOutput();
    }
    if (ftr.EnableMatrix()) {
        enableMatrixOutput();
    }
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

void Field::fillInitial() {
    t = 0;
    lastIterrationsCount = 0;
    if (transposed) {
        transpose();
    }

    for (size_t index = 0, len = width * height; index < len; ++index) {
        curr[index] = ftr.TStart();
    }
}

void Field::transpose() {
    transpose(transposed ? curr : prev);
    
    std::swap(hX, hY);
    std::swap(mySX, mySY);
    transposed = transposed == false;
}

void Field::nextTimeLayer() {
    std::swap(curr, prev);
    t += dT;
}

void Field::fillFactors(size_t row, bool first) {
    double *aF = maF + row * width;
    double *bF = mbF + row * width;
    double *cF = mcF + row * width;
    double *fF = mfF + row * width;

    double *rw = prev + row * width;
    double *brw = first ? rw : (curr + row * width);

    double lm0 = ftr.lambda(brw[0]), lmh = ftr.lambda(brw[1]);
    aF[0] = 0;
    cF[0] = dT * (lm0 + lmh) + hX * hX * ftr.ro(brw[0]) * ftr.cEf(brw[0]);
    bF[0] = -dT * (lm0 + lmh);
    fF[0] = hX * hX * ftr.ro(brw[0]) * ftr.cEf(brw[0]) * rw[0];

    double TPrev = brw[width - 1];
    double TPrev4 = TPrev * TPrev * TPrev * TPrev;
    double lmXX = ftr.lambda(brw[width - 1]), lmXXm1 = ftr.lambda(brw[width - 2]);
    aF[width - 1] = -dT * (lmXXm1 + lmXX);
    cF[width - 1] = dT * (lmXXm1 + lmXX)
                     + hX * hX * ftr.ro(brw[width - 1]) * ftr.cEf(brw[width - 1])
                     + 2 * hX * dT * ftr.alpha(t);
    bF[width - 1] = 0;
    fF[width - 1] = hX * hX * ftr.ro(brw[width - 1]) * ftr.cEf(brw[width - 1]) * rw[width - 1]
                     - 2 * hX * dT * ftr.sigma(t) * (TPrev4 - ftr.TEnv4())
                     + 2 * hX * dT * ftr.alpha(t) * ftr.TEnv();

    double lmXm1 = lm0, lmX = lmh, lmXp1;
    double mhh2rocdT;
    for (size_t index = 1; index < width - 1; ++index) {
        lmXp1 = ftr.lambda(brw[index + 1]);
        mhh2rocdT = - 2 * hX * hX * ftr.ro(brw[index]) * ftr.cEf(brw[index]) / dT;

        aF[index] = lmX + lmXm1;
        bF[index] = lmXp1 + lmX;
        cF[index] = -(lmXp1 + 2 * lmX + lmXm1) + mhh2rocdT;
        fF[index] = mhh2rocdT * rw[index];

        lmXm1 = lmX;
        lmX = lmXp1;
    }
}

void Field::firstPass(size_t row) {
    double *aF = maF + row * width;
    double *bF = mbF + row * width;
    double *cF = mcF + row * width;
    double *fF = mfF + row * width;

    double m = 0;
    for (size_t i = 1; i < width; ++i) {
        m = aF[i] / cF[i - 1];
        cF[i] -= m * bF[i - 1];
        fF[i] -= m * fF[i - 1];
    }
}

double Field::secondPass(size_t row, bool first) {
    double *bF = mbF + row * width;
    double *cF = mcF + row * width;
    double *fF = mfF + row * width;

    double *y = curr + row * width;
    double *py = first ? (prev + row * width) : y;

    double newValue = 0, maxDelta = 0;
    newValue = fF[width - 1] / cF[width - 1];
    maxDelta = fabs(newValue - py[width - 1]);
    y[width - 1] = newValue;

    for (ssize_t i = width - 2; i >= 0; --i) {
        newValue = (fF[i] - bF[i] * y[i + 1]) / cF[i];

        double newDelta = fabs(newValue - py[i]);
        maxDelta = std::max(maxDelta, newDelta);
        y[i] = newValue;
    }

    return maxDelta;
}

double Field::solve(size_t row, bool first) {
    firstPass(row);
    return secondPass(row, first);
}

void Field::test() {
    for (size_t index = 0; index < width * height; ++index) {
        prev[index] = myCoord * width * height + index;
    }

    MPI_Barrier(comm);
    transpose();
}

size_t Field::solveRows() {
    size_t maxIterationsCount = 0;

    for (size_t row = 0; row < height; ++row) {
        fillFactors(row, true);
        double delta = solve(row, true);
        size_t iterationsCount = 1;

        while (delta > epsilon) {
            fillFactors(row, false);
            delta = solve(row, false);
            ++iterationsCount;

            if (iterationsCount > MAX_ITTERATIONS_COUNT) {
                break;
            }
        }
        maxIterationsCount = std::max(maxIterationsCount, iterationsCount);
    }

    return maxIterationsCount;
}

void Field::solve() {
    if (done()) {
        return;
    }

    lastIterrationsCount = 0;

    nextTimeLayer();
    lastIterrationsCount += solveRows();

    nextTimeLayer();
    transpose();
    lastIterrationsCount += solveRows();
    transpose();

    printAll();
}

double Field::time() {
    return t;
}

bool Field::done() {
    return t >= ftr.TMax();
}

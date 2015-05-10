//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#include "field.h"
#include "factors.h"
#include <cmath>

int const MASTER = 0;
int const WAITER = 0;
int const NOBODY = MPI_PROC_NULL;
int const NOTHING = -1;

static size_t const MAX_ITTERATIONS_COUNT = 50;

Field::Field() {
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

    size_t maxDim = std::max(width, height);
    aF = new double[maxDim];
    bF = new double[maxDim];
    cF = new double[maxDim];
    fF = new double[maxDim];

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

    delete[] aF;
    delete[] bF;
    delete[] cF;
    delete[] fF;
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

void Field::transpose(double *arr) {
    for (size_t index = 0, len = width * height; index < len; ++index) {
        size_t newIndex = (index % width) * height + index / width;
        buff[newIndex] = arr[index];
    }
    for (size_t index = 0, len = width * height; index < len; ++index) {
        arr[index] = buff[index];
    }
}

void Field::transpose() {
    sendReceivePrevRows();
    transpose(prev);
    transpose(curr);
    
    std::swap(width, height);
    std::swap(hX, hY);
    transposed = transposed == false;

    std::swap(leftN, topN);
    std::swap(rightN, bottomN);
    std::swap(colComm, rowComm);
}

void Field::nextTimeLayer() {
    std::swap(curr, prev);
    t += dT;
}

void Field::fillFactors(size_t row, bool first) {
    double *rw = prev + row * width;
    double *brw = first && transposed == false ? rw : (curr + row * width);

    if (leftN == NOBODY) {
        double lm0 = ftr.lambda(brw[0]), lmh = ftr.lambda(brw[1]);
        aF[0] = 0;
        cF[0] = dT * (lm0 + lmh) + hX * hX * ftr.ro(brw[0]) * ftr.cEf(brw[0]);
        bF[0] = -dT * (lm0 + lmh);
        fF[0] = hX * hX * ftr.ro(brw[0]) * ftr.cEf(brw[0]) * rw[0];
    }

    if (rightN == NOBODY) {
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
    }

    for (size_t index = 1; index < width - 1; ++index) {
        double roc = ftr.ro(brw[index]) * ftr.cEf(brw[index]);
        double lmXm1 = ftr.lambda(brw[index - 1]), lmX = ftr.lambda(brw[index]), lmXp1 = ftr.lambda(brw[index + 1]);

        aF[index] = dT * (lmX + lmXm1);
        bF[index] = dT * (lmXp1 + lmX);
        cF[index] = -dT * (lmXp1 + 2 * lmX + lmXm1) - 2 * hX * hX * roc;
        fF[index] = -2 * hX * hX * roc * rw[index];
    }
}

double Field::solve(size_t row, bool first) {
    receiveFirstPass(row);
    double m = 0;
    for (size_t i = 1; i < width; ++i) {
        m = aF[i] / cF[i - 1];
        cF[i] -= m * bF[i - 1];
        fF[i] -= m * fF[i - 1];
    }
    sendFirstPass(row);

    double *y = curr + row * width;
    double *py = first ? (prev + row * width) : y;

    double newValue = 0, maxDelta = 0;
    if (rightN == NOBODY) {
        newValue = fF[width - 1] / cF[width - 1];
        maxDelta = fabs(newValue - py[width - 1]);
        y[width - 1] = newValue;
    } else {
        receiveSecondPass(row);
    }

    for (ssize_t i = width - 2; i >= 0; --i) {
        newValue = (fF[i] - bF[i] * y[i + 1]) / cF[i];

        double newDelta = fabs(newValue - py[i]);
        maxDelta = std::max(maxDelta, newDelta);
        y[i] = newValue;
    }
    sendSecondPass(row);

    reduceMaxDelta(maxDelta);
    sendReceiveCurrRowBorders(row);

    return maxDelta;
}

void Field::test() {
    size_t w = width;
    aF[0] = 0; aF[1] = 1; aF[2] = 1; aF[3] = 1; aF[4] = 1;
    bF[0] = 3; bF[1] = -1; bF[2] = 1; bF[3] = -1; bF[4] = 0;
    cF[0] = 4; cF[1] = 2; cF[2] = 4; cF[3] = 2; cF[4] = 2;
    fF[0] = 4; fF[1] = 2; fF[2] = 7.5; fF[3] = 1; fF[4] = -3;
    prev[0] = -0.540909; prev[1] = 2.05455; prev[2] = 1.56818; prev[3] = -0.827273; prev[4] = -1.08636;

    width = 5;

    solve(0, true);

    for (size_t i = 0; i < 5; ++i) {
        if (fabs(prev[i] - curr[i]) > 0.00001) {
            printf("Error: %lu\n", i);
        }
    }

    width = w;
}

size_t Field::solveRows() {
    size_t maxIterationsCount = 0;

    for (size_t row = (topN != NOBODY ? 1 : 0); row < height - (bottomN != NOBODY ? 1 : 0); ++row) {
        fillFactors(row, true);
        double delta = solve(row, true);
        size_t iterationsCount = 1;

        while (delta > epsilon) {
            fillFactors(row, false);
            delta = solve(row, false);
            ++iterationsCount;

            /*if (iterationsCount > MAX_ITTERATIONS_COUNT / 2) {
                printf("Warning! Iterations (%lu:%lu) on layer going to maximum\n",
                       transposed ? 0 : row, transposed ? row : 0);
            }*/

            if (iterationsCount > MAX_ITTERATIONS_COUNT) {
                /*printf("Error! Iterations (%lu:%lu) on layer [t: %.3f, index: %lu, delta: %.5f] achieve maximum\n",
                       transposed ? 0 : row, transposed ? row : 0, t, (unsigned long)(t / dT), delta);*/
                //exit(1);
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

    printMatrix();
    printViews();
    print();
}

double Field::time() {
    return t;
}

bool Field::done() {
    return t >= ftr.totalTime();
}

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

    size_t wh = width * height;
    prev = new double[wh];
    curr = new double[wh];
    buff = new double[wh];
    sendBuff = new double[wh];
    recvBuff = new double[wh];
    views = new double[ftr.ViewCount()];

    aF = new double[wh];
    bF = new double[wh];
    cF = new double[wh];
    fF = new double[wh];

    rowDeltas = new double[std::max(width, height)];

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
    delete[] sendBuff;
    delete[] recvBuff;
    delete[] views;

    delete[] aF;
    delete[] bF;
    delete[] cF;
    delete[] fF;

    delete[] rowDeltas;
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

    std::swap(leftN, topN);
    std::swap(rightN, bottomN);
    std::swap(colComm, rowComm);

    transposed = transposed == false;
}

void Field::nextTimeLayer() {
    std::swap(curr, prev);
    t += dT;
}

void Field::resetDeltas() {
    for (size_t i = 0; i < height; ++i) {
        rowDeltas[i] = __DBL_MAX__;
    }
}

void Field::fillFactors(bool first) {
    for (size_t row = (topN != NOBODY ? 1 : 0), len = height - (bottomN != NOBODY ? 1 : 0); row < len; ++row) {
        /*if (rowDeltas[row] <= epsilon) {
            continue;
        }*/

        double *raF = aF + row * width;
        double *rbF = bF + row * width;
        double *rcF = cF + row * width;
        double *rfF = fF + row * width;

        double *rw = prev + row * width;
        double *brw = first && transposed == false ? rw : (curr + row * width);

        if (leftN == NOBODY) {
            double lm0 = ftr.lambda(brw[0]), lmh = ftr.lambda(brw[1]);
            raF[0] = 0;
            rcF[0] = dT * (lm0 + lmh) + hX * hX * ftr.ro(brw[0]) * ftr.cEf(brw[0]);
            rbF[0] = -dT * (lm0 + lmh);
            rfF[0] = hX * hX * ftr.ro(brw[0]) * ftr.cEf(brw[0]) * rw[0];
        }

        if (rightN == NOBODY) {
            double TPrev = brw[width - 1];
            double TPrev4 = TPrev * TPrev * TPrev * TPrev;
            double lmXX = ftr.lambda(brw[width - 1]), lmXXm1 = ftr.lambda(brw[width - 2]);
            raF[width - 1] = -dT * (lmXXm1 + lmXX);
            rcF[width - 1] = dT * (lmXXm1 + lmXX)
                             + hX * hX * ftr.ro(brw[width - 1]) * ftr.cEf(brw[width - 1])
                             + 2 * hX * dT * ftr.alpha(t);
            rbF[width - 1] = 0;
            rfF[width - 1] = hX * hX * ftr.ro(brw[width - 1]) * ftr.cEf(brw[width - 1]) * rw[width - 1]
                             - 2 * hX * dT * ftr.sigma(t) * (TPrev4 - ftr.TEnv4())
                             + 2 * hX * dT * ftr.alpha(t) * ftr.TEnv();
        }

        for (size_t index = 1; index < width - 1; ++index) {
            double roc = ftr.ro(brw[index]) * ftr.cEf(brw[index]);
            double lmXm1 = ftr.lambda(brw[index - 1]), lmX = ftr.lambda(brw[index]), lmXp1 = ftr.lambda(brw[index + 1]);

            raF[index] = dT * (lmX + lmXm1);
            rbF[index] = dT * (lmXp1 + lmX);
            rcF[index] = -dT * (lmXp1 + 2 * lmX + lmXm1) - 2 * hX * hX * roc;
            rfF[index] = -2 * hX * hX * roc * rw[index];
        }
    }
}

void Field::firstPass() {
    receiveFirstPass();

    for (size_t row = (topN != NOBODY ? 1 : 0), len = height - (bottomN != NOBODY ? 1 : 0); row < len; ++row) {
        /*if (rowDeltas[row] <= epsilon) {
            continue;
        }*/

        double *raF = aF + row * width;
        double *rbF = bF + row * width;
        double *rcF = cF + row * width;
        double *rfF = fF + row * width;

        double m = 0;
        for (size_t i = 1; i < width; ++i) {
            m = raF[i] / rcF[i - 1];
            rcF[i] -= m * rbF[i - 1];
            rfF[i] -= m * rfF[i - 1];
        }
    }

    sendFirstPass();
}

void Field::secondPass(bool first) {
    receiveSecondPass();

    for (size_t row = (topN != NOBODY ? 1 : 0), len = height - (bottomN != NOBODY ? 1 : 0); row < len; ++row) {
        /*if (first == false && rowDeltas[row] <= epsilon) {
            continue;
        }*/

        double *rbF = bF + row * width;
        double *rcF = cF + row * width;
        double *rfF = fF + row * width;

        double *y = curr + row * width;
        double *py = first ? (prev + row * width) : y;

        double newValue = 0, maxDelta = 0;
        if (rightN == NOBODY) {
            newValue = rfF[width - 1] / rcF[width - 1];
            maxDelta = fabs(newValue - py[width - 1]);
            y[width - 1] = newValue;
        }

        for (ssize_t i = width - 2, minI = (leftN != NOBODY ? 1 : 0); i >= minI; --i) {
            newValue = (rfF[i] - rbF[i] * y[i + 1]) / rcF[i];

            double newDelta = fabs(newValue - py[i]);
            maxDelta = std::max(maxDelta, newDelta);
            y[i] = newValue;
        }
        
        rowDeltas[row] = maxDelta;
    }

    sendSecondPass();
    sendReceiveLeftBorders();
}

size_t Field::solveRows() {
    resetDeltas();

    fillFactors(true);
    firstPass();
    secondPass(true);
    double delta = reduceMaxDelta();
    size_t iterationsCount = 1;

    while (delta > epsilon) {
        fillFactors(false);
        firstPass();
        secondPass(false);
        delta = reduceMaxDelta();
        ++iterationsCount;

        if (iterationsCount > MAX_ITTERATIONS_COUNT) {
            /*printf("Error! Iterations (%lu:%lu) on layer [t: %.3f, index: %lu, delta: %.5f] achieve maximum\n",
             transposed ? 0 : row, transposed ? row : 0, t, (unsigned long)(t / dT), delta);*/
            //exit(1);
            break;
        }
    }

    return iterationsCount;
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

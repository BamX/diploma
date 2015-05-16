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
int const SEND_PACK_SIZE = 6;

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

    maF = new double[height * width];
    mbF = new double[height * width];
    mcF = new double[height * width];
    mfF = new double[height * width];
    calculatingRows = new bool[width];
    prevCalculatingRows = new bool[width];


    sendBuff = new double[width * SEND_PACK_SIZE];
    boolSendBuff = new bool[width];
    receiveBuff = new double[width * SEND_PACK_SIZE];

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
    delete[] calculatingRows;
    delete[] prevCalculatingRows;

    delete[] sendBuff;
    delete[] boolSendBuff;
    delete[] receiveBuff;
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
    transpose(transposed ? curr : prev);
    
    std::swap(hX, hY);
    std::swap(mySX, mySY);
    std::swap(topN, leftN);
    std::swap(bottomN, rightN);
    std::swap(width, height);
    transposed = transposed == false;
}

void Field::nextTimeLayer() {
    std::swap(curr, prev);
    t += dT;
}

void Field::resetCalculatingRows() {
    for (size_t index = 0; index < height; ++index) {
        calculatingRows[index] = prevCalculatingRows[index] = true;
    }
}

void Field::fillFactors(size_t row, bool first) {
    double *aF = maF + row * width;
    double *bF = mbF + row * width;
    double *cF = mcF + row * width;
    double *fF = mfF + row * width;

    double *rw = prev + row * width;
    double *brw = first ? rw : (curr + row * width);

    double lm0 = ftr.lambda(brw[0]), lmh = ftr.lambda(brw[1]);
    if (leftN == NOBODY) {
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

    double lmXm1 = lm0, lmX = lmh, lmXp1;
    for (size_t index = 1; index < width - 1; ++index) {
        lmXp1 = ftr.lambda(brw[index + 1]);
        double roc = ftr.ro(brw[index]) * ftr.cEf(brw[index]);

        aF[index] = dT * (lmX + lmXm1);
        bF[index] = dT * (lmXp1 + lmX);
        cF[index] = -dT * (lmXp1 + 2 * lmX + lmXm1) - 2 * hX * hX * roc;
        fF[index] = -2 * hX * hX * roc * rw[index];

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
    if (rightN == NOBODY) {
        newValue = fF[width - 1] / cF[width - 1];
        maxDelta = fabs(newValue - py[width - 1]);
        y[width - 1] = newValue;
    }

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

    if (transposed) {
        resetCalculatingRows();
        bool first = true;
        bool solving = true;
        while (solving) {
            solving = false;

            sendRecieveCalculatingRows();
            size_t fromRow = 0;
            while (fromRow < height) {
                size_t bundleSize = 0;

                recieveFirstPass(fromRow); // calculatingRows x [prevCalculatingRows] + (b + c + f) x [calculatingRows]
                size_t row = fromRow;
                for (; row < height; ++row) {
                    if (calculatingRows[row] == false) {
                        continue;
                    }
                    else {
                        solving = true;
                    }

                    fillFactors(row, first);
                    firstPass(row);

                    ++bundleSize;
                    if (bundleSize >= bundleSizeLimit) {
                        break;
                    }
                }
                sendFistPass(fromRow); // calculatingRows x [prevCalculatingRows] + (b + c + f) x [calculatingRows]
                fromRow = row;
            }

            if (solving == false) {
                break;
            }

            std::swap(prevCalculatingRows, calculatingRows);

            fromRow = 0;
            while (fromRow < height) {
                size_t bundleSize = 0;

                recieveSecondPass(fromRow); // (calculatingRows + y) x [prevCalculatingRows]
                size_t row = fromRow;
                for (; row < height; ++row) {
                    if (prevCalculatingRows[row] == false) {
                        calculatingRows[row] = false;
                        continue;
                    }

                    double delta = secondPass(row, first);

                    calculatingRows[row] = (rightN == NOBODY ? false : calculatingRows[row]) || delta > epsilon;

                    ++bundleSize;
                    if (bundleSize >= bundleSizeLimit) {
                        break;
                    }
                }
                sendSecondPass(fromRow); // (calculatingRows + y) x [prevCalculatingRows]

                fromRow = row;
            }

            first = false;
            ++maxIterationsCount;
            if (maxIterationsCount >= MAX_ITTERATIONS_COUNT) {
                solving = false;
                break;
            }
        }
    }
    else {
        size_t startRow = topN == NOBODY ? 0 : 1;
        size_t realHeight = height - (bottomN == NOBODY ? 0 : 1);
        for (size_t row = startRow; row < realHeight; ++row) {
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
    sendRecieveRows();
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

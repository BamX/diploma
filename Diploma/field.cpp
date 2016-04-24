//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#include "field.h"
#include "algo.h"
#include <cmath>

int const MASTER = 0;
int const WAITER = 3;
int const NOBODY = MPI_PROC_NULL;
int const NOTHING = -1;
int const SEND_PACK_SIZE = 6;

static size_t const MAX_ITTERATIONS_COUNT = 50;

Field::Field() {
    initFactors();

    width = algo::ftr().X1SplitCount();
    height = algo::ftr().X2SplitCount();

    hX = algo::ftr().X1() / (width - 1);
    hY = algo::ftr().X2() / (height - 1);
    dT = algo::ftr().totalTime() / algo::ftr().TimeSplitCount();
    epsilon = algo::ftr().Epsilon();
    transposed = false;
    fout = NULL;
    mfout = NULL;

    calculateNBS();

    prev = new double[height * width];
    curr = new double[height * width];
    buff = new double[height * width];
    views = new double[algo::ftr().ViewCount()];

    maF = new double[height * width];
    mbF = new double[height * width];
    mcF = new double[height * width];
    mfF = new double[height * width];
    calculatingRows = new bool[width];
    nextCalculatingRows = new bool[width];

    lastIterationsCount = lastWaitingCount = 0;

    sendBuff = new double[width * SEND_PACK_SIZE];
    boolSendBuff = new bool[width];
    receiveBuff = new double[width * SEND_PACK_SIZE];

    if (algo::ftr().EnablePlot()) {
        enablePlotOutput();
    }
    if (algo::ftr().EnableMatrix()) {
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
    delete[] nextCalculatingRows;

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
        curr[index] = algo::ftr().TStart();
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
        calculatingRows[index] = nextCalculatingRows[index] = true;
    }
}

void Field::fillFactors(size_t row, bool first) {
    double *aF = maF + row * width;
    double *bF = mbF + row * width;
    double *cF = mcF + row * width;
    double *fF = mfF + row * width;

    double *rw = prev + row * width;
    double *brw = first ? rw : (curr + row * width);

    algo::fillFactors(rw, brw, width, aF, bF, cF, fF, t, hX, dT, leftN == NOBODY, rightN == NOBODY);
}

void Field::firstPass(size_t row) {
    double *aF = maF + row * width;
    double *bF = mbF + row * width;
    double *cF = mcF + row * width;
    double *fF = mfF + row * width;

    algo::firstPass(width, aF, bF, cF, fF);
}

double Field::secondPass(size_t row, bool first) {
    double *bF = mbF + row * width;
    double *cF = mcF + row * width;
    double *fF = mfF + row * width;

    double *y = curr + row * width;
    double *py = first ? (prev + row * width) : y;

    double maxDelta = 0;
    algo::secondPass(py, y, width, bF, cF, fF, rightN == NOBODY, &maxDelta);

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

size_t Field::firstPasses(size_t fromRow, bool first, bool async) {
    if (async && checkIncomingFirstPass(fromRow) == false) {
        return 0;
    }

    recieveFirstPass(fromRow, first); // calculatingRows x [prevCalculatingRows] + (b + c + f) x [calculatingRows]
    size_t row = fromRow;
    for (size_t bundleSize = 0; row < height && bundleSize < bundleSizeLimit; ++row, ++bundleSize) {
        if (calculatingRows[row] == false) {
            continue;
        }

        firstPass(row);
    }
    sendFistPass(fromRow); // calculatingRows x [prevCalculatingRows] + (b + c + f) x [calculatingRows]

    return row;
}

size_t Field::secondPasses(size_t fromRow, bool first, bool async) {
    if (checkIncomingSecondPass(fromRow) == false) {
        if (async) {
            return 0;
        }
        else if (myCoord == 0) {
            ++lastWaitingCount;
        }
    }

    if (myCoord == 0) {
        ++lastIterationsCount;
    }

    recieveSecondPass(fromRow); // (nextCalculatingRows + y) x [prevCalculatingRows]

    size_t row = fromRow;
    for (size_t bundleSize = 0; row < height && bundleSize < bundleSizeLimit; ++row, ++bundleSize) {
        if (calculatingRows[row] == false) {
            nextCalculatingRows[row] = false;
            continue;
        }

        double delta = secondPass(row, first);

        nextCalculatingRows[row] = (rightN == NOBODY ? false : nextCalculatingRows[row]) || delta > epsilon;
    }
    sendSecondPass(fromRow); // (nextCalculatingRows + y) x [prevCalculatingRows]

    return row;
}

void Field::balanceBundleSize() {
    //printf("%zu\t%zu\n", lastWaitingCount, lastIterationsCount);
    if (lastWaitingCount > lastIterationsCount * algo::ftr().BalanceFactor()) {
        bundleSizeLimit = std::max(bundleSizeLimit - 1, algo::ftr().MinimumBundle());
    }
    else {
        bundleSizeLimit = std::min(bundleSizeLimit + 1, height / 2);
    }
    lastIterationsCount = 0;
    lastWaitingCount = 0;
}

size_t Field::solveRows() {
    size_t maxIterationsCount = 0;

    if (transposed) {
        resetCalculatingRows();
        bool first = true;

        if (myCoord == 0) {
            balanceBundleSize();
        }

        bool solving = true;
        while (solving) {
            sendRecieveCalculatingRows();
            if (first == false) {
                solving = false;
                for (size_t row = 0; row < height; ++row) {
                    if (calculatingRows[row]) {
                        solving = true;
                        break;
                    }
                }
                if (solving == false) {
                    break;
                }
            }

            for (size_t row = 0; row < height; ++row) {
                if (calculatingRows[row]) {
                    fillFactors(row, first);
                }
            }

            size_t fromFirstPassRow = 0;
            size_t fromSecondPassRow = 0;

            while (fromSecondPassRow < height) {
                size_t nextSecondPassRow = 0;
                if (fromSecondPassRow < height && fromFirstPassRow > fromSecondPassRow) {
                    nextSecondPassRow = secondPasses(fromSecondPassRow, first, true);
                }

                size_t nextFirstPassRow = 0;
                if (fromFirstPassRow < height) {
                    nextFirstPassRow = firstPasses(fromFirstPassRow, first, true);
                }

                if (nextFirstPassRow == 0 && nextSecondPassRow == 0) {
                    if (fromFirstPassRow < height) {
                        nextFirstPassRow = firstPasses(fromFirstPassRow, first, false);
                    }
                    else {
                        nextSecondPassRow = secondPasses(fromSecondPassRow, first, false);
                    }
                }

                if (nextFirstPassRow > 0) {
                    fromFirstPassRow = nextFirstPassRow;
                }
                if (nextSecondPassRow > 0) {
                    fromSecondPassRow = nextSecondPassRow;
                }
            }

            first = false;
            std::swap(nextCalculatingRows, calculatingRows);
            ++maxIterationsCount;
            if (maxIterationsCount >= MAX_ITTERATIONS_COUNT) {
                solving = false;
                break;
            }
        }
    }
    else {
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
    return t >= algo::ftr().TMax();
}

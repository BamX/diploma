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

size_t const MAX_ITTERATIONS_COUNT = 50;

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

size_t Field::solveRows() {
    return 0;
}

void Field::test() {
    for (size_t index = 0; index < width * height; ++index) {
        prev[index] = myCoord * width * height + index;
    }

    MPI_Barrier(comm);
    transpose();
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

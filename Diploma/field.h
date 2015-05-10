//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#ifndef __Diploma__field__
#define __Diploma__field__

#include "factors.h"
#include <fstream>
#include <mpi.h>

class Field {
    Factors ftr;
    std::ofstream *fout, *mfout;
    double nextFrameTime;

    double *prev, *curr, *buff, *views;
    double *aF, *bF, *cF, *fF;

    size_t width;
    size_t height;
    bool transposed;

    double t;
    double hX, hY, dT;
    double epsilon;
    size_t lastIterrationsCount;

    int myId;
    MPI_Comm comm;
    size_t mySX, mySY;
    int leftN, rightN, topN, bottomN;

    void calculateNBS();
    void calculateGrid(int numProcs, double stExpected, int &stX, int &stY);

    void fillFactors(size_t row, bool first);
    double solve(size_t row, bool first);
    size_t solveRows();
    void transpose(double *arr);
    void transpose();
    void nextTimeLayer();

    void enablePlotOutput();
    void enableMatrixOutput();

    void print();
    void printMatrix();
    void printViews();
    void debug(const char *name);

    void sendReceivePrevRows();
    void sendReceiveCurrRowBorders(size_t row);
    void sendFirstPass(size_t row);
    void receiveFirstPass(size_t row);
    void sendSecondPass(size_t row);
    void receiveSecondPass(size_t row);
    void reduceMaxDelta(double &maxDelta);
    void reduceViews();

public:
    Field();
    ~Field();

    void randomFill();
    void test();

    void fillInitial();
    void solve();
    double time();
    bool done();

    double view(double x1, double x2);
    double view(size_t index);

};

#endif /* defined(__Diploma__field__) */

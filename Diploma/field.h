//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#ifndef __Diploma__field__
#define __Diploma__field__

#include "factors.h"
#include <fstream>
#include <mpi.h>

extern int const MASTER;
extern int const WAITER;
extern int const NOBODY;
extern int const NOTHING;

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

    int myId, numProcs, myCoord;
    MPI_Comm comm;
    MPI_Datatype mpiAllType;
    size_t mySX, mySY;
    int topN, bottomN;

    void calculateNBS();

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

    void reduceViews();

public:
    Field();
    ~Field();
    void finalize();

    void test();
    void testPrint();

    void fillInitial();
    void solve();
    double time();
    bool done();

    double view(double x1, double x2);
    double view(size_t index);

};

#endif /* defined(__Diploma__field__) */

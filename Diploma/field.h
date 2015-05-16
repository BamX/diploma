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
extern int const SEND_PACK_SIZE;

class Field {
    Factors ftr;
    std::ofstream *fout, *mfout;
    double nextFrameTime;

    double *prev, *curr, *buff, *views;
    double *maF, *mbF, *mcF, *mfF;
    bool *calculatingRows, *prevCalculatingRows;

    size_t width;
    size_t height;
    bool transposed;

    double t;
    double hX, hY, dT;
    double epsilon;
    size_t lastIterrationsCount;

    int myId, numProcs, myCoord;
    MPI_Comm comm;
    size_t mySX, mySY;
    int topN, bottomN, leftN, rightN;
    size_t bundleSizeLimit;
    double *sendBuff, *receiveBuff;
    bool *boolSendBuff;

    void calculateNBS();

    void fillFactors(size_t row, bool first);
    void firstPass(size_t row);
    double secondPass(size_t row, bool first);
    double solve(size_t row, bool first);

    size_t solveRows();
    void transpose(double *arr);
    void transpose();
    void nextTimeLayer();
    void resetCalculatingRows();

    void enablePlotOutput();
    void enableMatrixOutput();

    void print();
    void printMatrix();
    void printViews();
    void debug(const char *name);

    void sendRecieveRows();
    void sendRecieveCalculatingRows();
    void sendFistPass(size_t fromRow);
    void recieveFirstPass(size_t fromRow);
    void sendSecondPass(size_t fromRow);
    void recieveSecondPass(size_t fromRow);

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

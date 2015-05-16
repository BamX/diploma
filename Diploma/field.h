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
    float nextFrameTime;

    float *prev, *curr, *buff, *views;
    float *maF, *mbF, *mcF, *mfF;
    bool *calculatingRows, *prevCalculatingRows;

    size_t width;
    size_t height;
    bool transposed;

    float t;
    float hX, hY, dT;
    float epsilon;
    size_t lastIterrationsCount;

    int myId, numProcs, myCoord;
    MPI_Comm comm;
    size_t mySX, mySY;
    int topN, bottomN, leftN, rightN;
    size_t bundleSizeLimit;
    float *sendBuff, *receiveBuff;
    bool *boolSendBuff;

    void calculateNBS();

    void fillFactors(size_t row, bool first);
    void firstPass(size_t row);
    float secondPass(size_t row, bool first);
    float solve(size_t row, bool first);

    size_t solveRows();
    void transpose(float *arr);
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
    float time();
    bool done();

    float view(float x1, float x2);
    float view(size_t index);

};

#endif /* defined(__Diploma__field__) */

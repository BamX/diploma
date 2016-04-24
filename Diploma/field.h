//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#ifndef __Diploma__field__
#define __Diploma__field__

#include <fstream>
#include <mpi.h>

extern int const MASTER;
extern int const WAITER;
extern int const NOBODY;
extern int const NOTHING;
extern int const SEND_PACK_SIZE;

class Field {
protected:
    std::ofstream *fout, *mfout;
    double nextFrameTime;

    double *prev, *curr, *buff, *views;
    double *maF, *mbF, *mcF, *mfF;
    bool *calculatingRows, *nextCalculatingRows;

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
    size_t lastWaitingCount, lastIterationsCount;

    void initFactors();
    void calculateNBS();

    void fillFactors(size_t row, bool first);
    void firstPass(size_t row);
    double secondPass(size_t row, bool first);
    size_t firstPasses(size_t fromRow, bool first, bool async);
    size_t secondPasses(size_t fromRow, bool first, bool async);
    double solve(size_t row, bool first);

    size_t solveRows();
    void transpose(double *arr);
    void transpose();
    void nextTimeLayer();
    void resetCalculatingRows();

    void enablePlotOutput();
    void enableMatrixOutput();

    void printAll();
    void printConsole();
    void printMatrix();
    void printViews();
    void debug(const char *name);

    void sendRecieveCalculatingRows();
    void balanceBundleSize();

    void sendFistPass(size_t fromRow);
    bool checkIncomingFirstPass(size_t fromRow);
    void recieveFirstPass(size_t fromRow, bool first);
    void sendSecondPass(size_t fromRow);
    bool checkIncomingSecondPass(size_t fromRow);
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

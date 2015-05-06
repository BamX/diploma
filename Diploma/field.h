//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#ifndef __Diploma__field__
#define __Diploma__field__

#include "factors.h"
#include <fstream>

class Field {
    Factors ftr;
    std::ofstream *fout, *mfout;
    double nextFrameTime;

    double *prev, *curr, *buff;
    double *aF, *bF, *cF, *fF;

    size_t width;
    size_t height;
    bool transposed;

    double t;
    double hX, hY, dT;
    double epsilon;
    size_t lastIterrationsCount;

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

//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#ifndef __Diploma__field__
#define __Diploma__field__

#include <iostream>

class Field {
    double *data, *buff;
    double *aF, *bF, *cF, *fF;

    size_t width;
    size_t height;
    bool transposed;

    double t;
    double hX, hY, dT;
    double epsilon;

    void fillFactors(size_t row, bool first);
    double solve(size_t row, bool first);
    void solveRows();
    void transpose();
    void flushBuffer();

public:
    Field(size_t _width, size_t _height, size_t _tStep, double _epsilon = 0.0001);
    ~Field();

    inline double& at(size_t row, size_t col);

    void print();
    void randomFill();

    void fillInitial();
    void solve();
    double time();
    bool done();

};

#endif /* defined(__Diploma__field__) */

//
//  Copyright (c) 2015 BX23. All rights reserved.
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
    double epsilon;

    void fillFactors(size_t line, bool first);
    double solve(size_t row, bool first);
    void solveRows();
    void transpose();
    void flushBuffer();

public:
    Field(size_t _width, size_t _height, double _epsilon = 0.0001) {
        width = _width;
        height = _height;
        epsilon = _epsilon;
        transposed = false;

        data = new double[height * width];
        buff = new double[height * width];

        size_t maxDim = std::max(width, height);
        aF = new double[maxDim];
        bF = new double[maxDim];
        cF = new double[maxDim];
        fF = new double[maxDim];
    }

    ~Field() {
        delete[] data;
        delete[] buff;

        delete[] aF;
        delete[] bF;
        delete[] cF;
        delete[] fF;
    }

    inline double& at(size_t row, size_t col);

    void print();
    void randomFill();

    void fill(double val);
    void solve();

};

#endif /* defined(__Diploma__field__) */

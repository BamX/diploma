//
//  Copyright (c) 2015 BX23. All rights reserved.
//

#ifndef __Diploma__field__
#define __Diploma__field__

#include <stdio.h>

class Field {
    double *data;
    size_t width;
    size_t height;

public:
    Field(size_t width, size_t height) {
        this->width = width;
        this->height = height;

        this->data = new double[height * width];
    }

    ~Field() {
        delete[] data;
    }

    inline double& at(size_t row, size_t col);
    void print();

    void fill(double val);
};

#endif /* defined(__Diploma__field__) */

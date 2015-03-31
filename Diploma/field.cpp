//
//  Copyright (c) 2015 BX23. All rights reserved.
//

#include "field.h"

inline double& Field::at(size_t row, size_t col) {
    return data[row * height + col];
}

void Field::print() {
    for (size_t index = 0, len = width * height; index < len; ++index) {
        printf("%.2f\t", data[index]);
        if ((index + 1) % width == 0) {
            printf("\n");
        }
    }
}

void Field::fill(double val) {
    for (size_t index = 0, len = width * height; index < len; ++index) {
        data[index] = val;
    }
}

//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#include "field.h"
#include "factors.h"
#include <cmath>

Field::Field(size_t _width, size_t _height, size_t _tStep, double _epsilon) {
    width = _width;
    height = _height;

    hX = ftr::X1 / _width;
    hY = ftr::X2 / _height;
    dT = _tStep;
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

Field::~Field() {
    delete[] data;
    delete[] buff;

    delete[] aF;
    delete[] bF;
    delete[] cF;
    delete[] fF;
}

inline double& Field::at(size_t row, size_t col) {
    return data[row * width + col];
}

void Field::randomFill() {
    for (size_t index = 0, len = width * height; index < len; ++index) {
        buff[index] = data[index] = rand() % 100;
    }
}

void Field::print() {
    printf("Field [%zux%zu]:\n", width, height);
    for (size_t index = 0, len = width * height; index < len; ++index) {
        printf("%.2f\t", data[index]);
        if ((index + 1) % width == 0) {
            printf("\n");
        }
    }
}

void Field::fillInitial() {
    t = 0;
    for (size_t index = 0, len = width * height; index < len; ++index) {
        buff[index] = data[index] = ftr::TStart;
    }
}

void Field::transpose() {
    for (size_t index = 0, len = width * height; index < len; ++index) {
        size_t newIndex = (index % width) * height + index / width;
        buff[newIndex] = data[index];
    }
    for (size_t index = 0, len = width * height; index < len; ++index) {
        data[index] = buff[index];
    }
    std::swap(width, height);
    transposed = transposed == false;
}

void Field::flushBuffer() {
    for (size_t index = 0, len = width * height; index < len; ++index) {
        data[index] = buff[index];
    }
}

void Field::fillFactors(size_t line, bool first) {
    // TODO: fill factors
}

double Field::solve(size_t row, bool first)
{
    fillFactors(row, first);

    double m = 0;
    for (size_t i = 1; i < width; ++i) {
        m = aF[i] / cF[i - 1];
        cF[i] -= m * bF[i - 1];
        fF[i] -= m * fF[i - 1];
    }

    size_t indexPrefix = row * width;
    double newValue = fF[width - 1] / cF[width - 1];
    double maxDelta = fabs(newValue - data[indexPrefix + width - 1]);
    data[indexPrefix + width - 1] = newValue;

    for (ssize_t i = width - 2; i >= 0; --i) {
        newValue = (fF[i] - bF[i] * data[indexPrefix + i + 1]) / cF[i];
        maxDelta = std::max(maxDelta, fabs(newValue - data[indexPrefix + i]));
        data[indexPrefix + i] = newValue;
    }

    return maxDelta;
}

void Field::solveRows() {
    for (size_t row = 0; row < height; ++row) {
        double delta = solve(row, true);

        while (delta > epsilon) {
            delta = solve(row, false);
        }
    }
    flushBuffer();
}

void Field::solve() {
    if (done()) {
        return;
    }

    solveRows();
    transpose();
    solveRows();
    transpose();

    t += dT;
}

double Field::time() {
    return t;
}

bool Field::done() {
    return t >= ftr::totalTime;
}

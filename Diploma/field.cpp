//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#include "field.h"
#include "factors.h"
#include <iostream>
#include <cmath>

Field::Field(size_t _width, size_t _height, size_t _tLength, double _epsilon) {
    width = _width;
    height = _height;

    hX = ftr.X1() / _width;
    hY = ftr.X2() / _height;
    dT = ftr.totalTime() / _tLength;
    epsilon = _epsilon;
    transposed = false;
    fout = NULL;

    data = new double[height * width];
    buff = new double[height * width];

    size_t maxDim = std::max(width, height);
    aF = new double[maxDim];
    bF = new double[maxDim];
    cF = new double[maxDim];
    fF = new double[maxDim];
}

Field::~Field() {
    if (fout != NULL) {
        fout->close();
        delete fout;
    }
    
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

double Field::view() {
    size_t x1index = ftr.X1View() / hX;
    size_t x2index = ftr.X2View() / hY;
    double x1factor = ftr.X1View() - x1index * hX;
    double x2factor = ftr.X2View() - x2index * hY;

    double value = data[x1index * width + x2index] + x1factor * data[x1index * width + x2index + 1];
    value += x2factor * (data[(x1index + 1) * width + x2index] + x1factor * data[(x1index + 1) * width + x2index + 1]);

    return value;
}

void Field::enableFileOutput() {
    if (fout != NULL) {
        fout->close();
    }
    fout = new std::ofstream("view.csv");
}

void Field::print() {
    printf("Field [%zux%zu](itrs: %zu, time: %.5f)\tview: %.7f\n", width, height, lastIterrationsCount, t, view());
}

void Field::fillInitial() {
    t = 0;
    lastIterrationsCount = 0;
    if (transposed) {
        transpose();
    }

    for (size_t index = 0, len = width * height; index < len; ++index) {
        buff[index] = data[index] = ftr.TStart();
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
    std::swap(hX, hY);
    transposed = transposed == false;
}

void Field::flushBuffer() {
    std::swap(data, buff);
}

void Field::fillFactors(size_t row, bool first) {
    double *rw = data + row * width;
    double *brw = first ? rw : (buff + row * width);

    double TPrev = brw[width - 1];
    double TPrev4 = TPrev * TPrev * TPrev * TPrev;

    aF[0] = 0;
    cF[0] = 1;
    bF[0] = 1 / (hX + 1);
    fF[0] = 0;

    double lmXX = ftr.lambda(brw[width - 1]), lmXXm1 = ftr.lambda(brw[width - 2]);
    aF[width - 1] = -dT * (lmXXm1 + lmXX);
    cF[width - 1] = -dT * (lmXXm1 + lmXX) - hX * hX + 2 * hX * dT * ftr.alpha(t);
    bF[width - 1] = 0;
    fF[width - 1] = hX * hX * rw[width - 1]
                     - 2 * hX * dT * ftr.sigma(t) * (TPrev4 - ftr.TEnv4())
                     + 2 * hX * dT * ftr.alpha(t) * ftr.TEnv();

    for (size_t index = 1; index < width - 1; ++index) {
        double roc = ftr.ro(brw[index]) * ftr.cEf(brw[index]);
        double lmXm1 = ftr.lambda(brw[index - 1]), lmX = ftr.lambda(brw[index]), lmXp1 = ftr.lambda(brw[index + 1]);

        aF[index] = dT * (lmX + lmXm1);
        bF[index] = dT * (lmXp1 + lmX);
        cF[index] = dT * (lmXp1 + 2 * lmX + lmXm1) - 2 * hX * hX * roc;
        fF[index] = -2 * hX * hX * roc * rw[index];
    }
}

double Field::solve(size_t row)
{
    double m = 0;
    for (size_t i = 1; i < width; ++i) {
        m = aF[i] / cF[i - 1];
        cF[i] -= m * bF[i - 1];
        fF[i] -= m * fF[i - 1];
    }

    double *y = buff + row * width;

    double newValue = fF[width - 1] / cF[width - 1];
    double maxDelta = fabs(newValue - y[width - 1]);
    y[width - 1] = newValue;

    for (ssize_t i = width - 2; i >= 0; --i) {
        newValue = (fF[i] - bF[i] * y[i + 1]) / cF[i];

        maxDelta = std::max(maxDelta, fabs(newValue - y[i]));
        y[i] = newValue;
    }

    return maxDelta;
}

void Field::test() {
    size_t w = width;
    aF[0] = 0; aF[1] = 1; aF[2] = 1; aF[3] = 1; aF[4] = 1;
    bF[0] = 3; bF[1] = -1; bF[2] = 1; bF[3] = -1; bF[4] = 0;
    cF[0] = 4; cF[1] = 2; cF[2] = 4; cF[3] = 2; cF[4] = 2;
    fF[0] = 4; fF[1] = 2; fF[2] = 7.5; fF[3] = 1; fF[4] = -3;
    data[0] = -0.540909; data[1] = 2.05455; data[2] = 1.56818; data[3] = -0.827273; data[4] = -1.08636;

    width = 5;

    solve(0);

    for (size_t i = 0; i < 5; ++i) {
        if (fabs(data[i] - buff[i]) > 0.00001) {
            printf("Error: %lu\n", i);
        }
    }

    width = w;
}

size_t Field::solveRows() {
    size_t maxIterationsCount = 0;

    for (size_t row = 0; row < height; ++row) {
        fillFactors(row, true);
        double delta = solve(row);
        size_t iterationsCount = 1;

        while (delta > epsilon) {
            fillFactors(row, false);
            delta = solve(row);
            ++iterationsCount;
        }
        maxIterationsCount = std::max(maxIterationsCount, iterationsCount);
    }

    flushBuffer();
    return maxIterationsCount;
}

void Field::solve() {
    if (done()) {
        return;
    }

    lastIterrationsCount = 0;

    lastIterrationsCount += solveRows();
    transpose();

    lastIterrationsCount += solveRows();
    transpose();

    t += dT;

    if (fout != NULL) {
        *fout << t << "," << view() << "\n";
    }
}

double Field::time() {
    return t;
}

bool Field::done() {
    return t >= ftr.totalTime();
}

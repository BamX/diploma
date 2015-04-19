//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#include "field.h"
#include "factors.h"
#include <iostream>
#include <cmath>

Field::Field(size_t _width, size_t _height, size_t _tStep, double _epsilon) {
    width = _width;
    height = _height;

    hX = ftr.X1() / _width;
    hY = ftr.X2() / _height;
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
    size_t x1index = ftr.X1View() / hX;
    size_t x2index = ftr.X2View() / hY;
    double x1factor = ftr.X1View() - x1index * hX;
    double x2factor = ftr.X2View() - x2index * hY;

    double value = data[x1index * width + x2index] + x1factor * data[x1index * width + x2index + 1];
    value += x2factor * (data[(x1index + 1) * width + x2index] + x1factor * data[(x1index + 1) * width + x2index + 1]);

    printf("Field [%zux%zu](itrs: %zu, time: %.1f)\tview: %.2f\n", width, height, lastIterrationsCount, t, value);
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
    transposed = transposed == false;
}

void Field::flushBuffer() {
    for (size_t index = 0, len = width * height; index < len; ++index) {
        data[index] = buff[index];
    }
}

double Field::lambda(size_t row, size_t x) {
    return ftr.lambda(at(row, x));
}

double Field::a(size_t row, size_t x) {
    return 0.5 * (lambda(row, x) + lambda(row, x + 1));
}

double Field::roc(size_t row, size_t x) {
    return ftr.ro(at(row, x)) * ftr.cEf(at(row, x));
}

void Field::fillFactors(size_t row, bool first) {
    size_t indexPrefix = row * width;
    double h = transposed ? hY : hX;
    double thF = dT / (h * h);

    aF[0] = 0;
    cF[0] = 1;
    bF[0] = -(dT * lambda(row, 1) / (dT * lambda(row, 1) + h * h / 2));
    fF[0] = h * h / 2 * data[indexPrefix] / (dT * lambda(row, 1) + h * h / 2);

    double TPrev = data[indexPrefix + width - 1];
    double TPrev4 = TPrev * TPrev * TPrev * TPrev;
    aF[width - 1] = dT * lambda(row, width - 2) / (dT * ftr.alpha(t) * h + dT * lambda(row, width - 2) - h * h / 2);
    cF[width - 1] = 1;
    bF[width - 1] = 0;
    fF[width - 1] =
        dT * h * (ftr.alpha(t) * ftr.TEnv4() + h / (2 * dT) * data[indexPrefix + width - 1] - ftr.sigma(t) * TPrev4) /
        (dT * ftr.alpha(t) * h + dT * lambda(row, width - 2) - h * h / 2);

    for (size_t index = 1; index < width - 1; ++index) {
        double thFROC = thF / roc(row, index);
        fF[index] = -data[indexPrefix + index];
        aF[index] = thFROC * a(row, index + 1);
        bF[index] = thFROC * a(row, index - 1);
        cF[index] = -(1 + thFROC * (a(row, index + 1) - a(row, index)));
    }
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
    double maxDelta = fabs(newValue - buff[indexPrefix + width - 1]);
    buff[indexPrefix + width - 1] = newValue;

    for (ssize_t i = width - 2; i >= 0; --i) {
        newValue = (fF[i] - bF[i] * data[indexPrefix + i + 1]) / cF[i];
        maxDelta = std::max(maxDelta, fabs(newValue - buff[indexPrefix + i]));
        buff[indexPrefix + i] = newValue;
    }

    return maxDelta;
}

size_t Field::solveRows() {
    size_t maxIterationsCount = 0;

    for (size_t row = 0; row < height; ++row) {
        double delta = solve(row, true);
        size_t iterationsCount = 1;

        while (delta > epsilon) {
            delta = solve(row, false);
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
}

double Field::time() {
    return t;
}

bool Field::done() {
    return t >= ftr.totalTime();
}

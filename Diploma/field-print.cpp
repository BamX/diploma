//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#include "field.h"
#include <cmath>

double Field::view(double x1, double x2) {
    ssize_t x1index = floor(x1 / hX) - mySX + (leftN != NOBODY ? 1 : 0);
    ssize_t x2index = floor(x2 / hY) - mySY + (topN != NOBODY ? 1 : 0);

    bool notInMyX1 = x1index < (leftN != NOBODY ? 1 : 0) || x1index >= width - (rightN != NOBODY ? 1 : 0);
    bool notInMyX2 = x2index < (topN != NOBODY ? 1 : 0) || x2index >= height - (bottomN != NOBODY ? 1 : 0);
    if (notInMyX1 || notInMyX2) {
        return NOTHING;
    }

    /*double x1factor = x1 - x1index * hX;
     double x2factor = x2 - x2index * hY;

     double value = curr[x2index * width + x1index] + x1factor * curr[x2index * width + x1index + 1];
     value += x2factor * (curr[(x2index + 1) * width + x1index] + x1factor * curr[(x2index + 1) * width + x1index + 1]);*/

    return curr[x2index * width + x1index];
}

double Field::view(size_t index) {
    return view(ftr.X1View(index), ftr.X2View(index));
}

void Field::enablePlotOutput() {
    if (myId != MASTER) {
        return;
    }

    if (fout != NULL) {
        fout->close();
    }
    fout = new std::ofstream("view.csv");
}

void Field::enableMatrixOutput() {
    if (myId != MASTER) {
        return;
    }

    if (mfout != NULL) {
        mfout->close();
    }
    mfout = new std::ofstream("matrix.csv");
}

void Field::printMatrix() {
    // TODO: print without overlapses
    if (mfout != NULL && t > nextFrameTime) {
        nextFrameTime += ftr.totalTime() / ftr.MatrixFramesCount();

        size_t index = 0;
        for (size_t row = 0; row < height; ++row) {
            for (size_t col = 0; col < width; ++col, ++index) {
                *mfout << curr[index];
                if (col < width - 1) {
                    *mfout << " ";
                }
            }
            *mfout << "\n";
        }
        *mfout << "\n";
        mfout->flush();
    }
}

void Field::printViews() {
    if (ftr.EnablePlot()) {
        reduceViews();
    }
    if (fout != NULL) {
        *fout << t;
        for (size_t index = 0, len = ftr.ViewCount(); index < len; ++index) {
            *fout << "," << views[index];
        }
        *fout << "\n";
        fout->flush();
    }
}

void Field::print() {
    double viewValue = view(ftr.DebugView());

    if (fabs(viewValue - NOTHING) > __DBL_EPSILON__) {
        printf("Field[%d] (itrs: %zu, time: %.5f)\tview: %.7f\n",
               myId, lastIterrationsCount, t, viewValue);
    }
}

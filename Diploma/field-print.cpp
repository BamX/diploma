//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#include "field.h"
#include <cmath>

#include <sys/types.h>
#include <unistd.h>

double Field::view(double x1, double x2) {
    ssize_t x1index = floor(x1 / hX) - mySX;
    ssize_t x2index = floor(x2 / hY) - mySY;

    bool notInMyX1 = x1index < 0 || x1index >= width;
    bool notInMyX2 = x2index < 0 || x2index >= height;
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

void Field::testPrint() {
    MPI_Barrier(comm);
    sleep(1);

    for (int p = 0; p < numProcs; ++p) {
        if (p == myCoord) {
            for (int i = 0, index = 0; i < height; ++i) {
                for (int j = 0; j < width; ++j, ++index) {
                    printf("%.0f\t", prev[index]);
                }
                printf("\n");
            }
        }
        MPI_Barrier(comm);
        sleep(1);
    }
}

//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#include "field.h"
#include <cmath>

#include <sys/types.h>
#include <unistd.h>

double Field::view(double x1, double x2) {
    ssize_t x1index = floor(x1 / hX) - mySX;
    ssize_t x2index = floor(x2 / hY) - mySY + (topN != NOBODY ? 1 : 0);

    bool notInMyX1 = x1index < 0 || x1index >= width;
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
    if (mfout != NULL) {
        mfout->close();
    }
    mfout = new std::ofstream("matrix.csv", std::ios::trunc);
}

void Field::printAll() {
    if (t > nextFrameTime) {
        nextFrameTime += ftr.TMax() / ftr.FramesCount();

        printConsole();
        printViews();
        printMatrix();
    }
}

void Field::printConsole() {
    if (ftr.EnableConsole()) {
        double viewValue = view(ftr.DebugView());

        if (fabs(viewValue - NOTHING) > __DBL_EPSILON__) {
            printf("Field[%d] (itrs: %zu, bsL %zu, time: %.5f)\tview: %.7f\n",
                   myId, lastIterrationsCount, bundleSizeLimit, t, viewValue);
        }
    }
}

void Field::printViews() {
    if (ftr.EnablePlot()) {
        reduceViews();

        if (fout != NULL) {
            *fout << t;
            for (size_t index = 0, len = ftr.ViewCount(); index < len; ++index) {
                *fout << "," << views[index];
            }
            *fout << "\n";
            fout->flush();
        }
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

void Field::logState(const char *name) {
    if (ftr.EnableStates()) {
        auto elapsed = std::chrono::high_resolution_clock::now() - startTime;
        unsigned long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
        char buff[100];
        snprintf(buff, 100, "[%d](%llu):%s\n", myId, microseconds, name);
        stateMessages.push_back(buff);
    }
}

void Field::printStates() {
    if (ftr.EnableStates()) {
        for (int p = 0; p < numProcs; ++p) {
            if (p == myCoord) {
                std::ofstream *fout = NULL;
                if (p == 0) {
                    fout = new std::ofstream("states.txt");
                } else {
                    fout = new std::ofstream("states.txt", std::ios::app);
                }

                *fout <<
                    "[" << myId <<
                    "]:w " << width <<
                    ":h " << height <<
                    ":t " << topN <<
                    ":b " << bottomN << "\n";

                for (size_t index = 0; index < stateMessages.size(); ++index) {
                    *fout << stateMessages[index];
                }
                fout->close();
                delete fout;
                stateMessages.resize(0);
            }
            MPI_Barrier(comm);
        }
    }
}

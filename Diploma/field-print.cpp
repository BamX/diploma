//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#include "field.h"
#include "algo.h"
#include <cmath>

#include <sys/types.h>
#include <unistd.h>

double Field::view(double x1, double x2) {
    long x1index = floor(x1 / hX) - mySX;
    long x2index = floor(x2 / hY) - mySY + (topN != NOBODY ? 1 : 0);

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
    return view(algo::ftr().X1View(index), algo::ftr().X2View(index));
}

void Field::enablePlotOutput() {
    if (myId != MASTER) {
        return;
    }

    if (fout != NULL) {
        fout->close();
        delete fout;
    }
    fout = new std::ofstream(algo::ftr().PlotFilename());
}

void Field::enableMatrixOutput() {
    if (mfout != NULL) {
        mfout->close();
        delete mfout;
    }
    mfout = new std::ofstream("matrix.csv", std::ios::trunc);
}

void Field::enableBucketsOutput() {
    if (isBucketsMaster() == false) {
        return;
    }
    
    if (bfout != NULL) {
        bfout->close();
        delete bfout;
    }
    bfout = new std::ofstream(algo::ftr().BucketsFilename());
    for (size_t i = 0; i < numProcs; ++i) {
        *bfout << "n" << i;
        if (i < numProcs - 1) {
             *bfout << ",";
        }
    }
    *bfout << "\n";
}

void Field::enableWeightsOutput() {
    if (isBucketsMaster() == false) {
        return;
    }

    if (wfout != NULL) {
        wfout->close();
        delete wfout;
    }
    wfout = new std::ofstream(algo::ftr().WeightsFilename());
    for (size_t i = 0, len = weightsSize(); i < len; ++i) {
        *wfout << "r" << i;
        if (i < len - 1) {
            *wfout << ",";
        }
    }
    *wfout << "\n";
}

void Field::enableTimesOutput() {
    if (tfout != NULL) {
        tfout->close();
        delete tfout;
    }
    char nameBuff[255] = {0};
    sprintf(nameBuff, "%s.%lu.csv", algo::ftr().TimesFilenamePrefix().data(), (long)myCoord);
    tfout = new std::ofstream(nameBuff);
}

void Field::printAll() {
    if (t > nextFrameTime) {
        nextFrameTime += algo::ftr().TMax() / algo::ftr().FramesCount();

        printConsole();
        printViews();
        printMatrix();
    }
}

void Field::printConsole() {
}

void Field::printViews() {
    if (algo::ftr().EnablePlot()) {
        reduceViews();

        if (fout != NULL) {
            *fout << t;
            for (size_t index = 0, len = algo::ftr().ViewCount(); index < len; ++index) {
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

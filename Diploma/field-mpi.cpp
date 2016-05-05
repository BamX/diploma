//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#include "field.h"
#include "algo.h"
#include <cmath>
#include <sys/types.h>
#include <unistd.h>

#define DEBUG_PRINT
//#define DEBUG_WAIT

void Field::debug(const char *name) {
#ifdef DEBUG_PRINT
    printf("I'm %d before %s\n", myId, name);
#endif
}

void Field::initFactors() {
    int procs, id;
    MPI_Comm_size(MPI_COMM_WORLD, &procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    for (size_t p = 0; p < procs; ++p) {
        if (p == id) {
            algo::ftr().initFactors();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void Field::calculateNBS() {
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    int dims[] = { numProcs };
    int wrap[] = { 0 };
    MPI_Cart_create(MPI_COMM_WORLD, 1, dims, wrap, 1, &comm);

    MPI_Comm_rank(comm, &myId);
    MPI_Cart_coords(comm, myId, 1, &myCoord);
    MPI_Cart_shift(comm, 0, 1, &topN, &bottomN);
    rightN = leftN = NOBODY;

    size_t newHeight = height / numProcs;

    mySY = newHeight * myCoord;
    mySX = 0;

    if (bottomN == NOBODY) {
        newHeight = height - newHeight * (numProcs - 1);
    }

    height = newHeight;

#ifdef DEBUG_WAIT
    printf("I'm %d(%d)\n", myId, ::getpid());
    int waiter = myId;
    while (waiter == WAITER) sleep(5);
#endif
}

void Field::finalize() {
}

void Field::reduceViews() {
    for (size_t index = 0, len = algo::ftr().ViewCount(); index < len; ++index) {
        double value = view(index);
        double result = 0;
        MPI_Reduce(&value, &result, 1, MPI_DOUBLE, MPI_MAX, MASTER, comm);
        views[index] = result;
    }
}

void Field::printMatrix() {
    if (algo::ftr().EnableMatrix()) {
        for (int p = 0; p < numProcs; ++p) {
            if (p == myCoord) {
                if (mfout != NULL) {
                    mfout->close();
                }
                mfout = new std::ofstream("matrix.csv", std::ios::app);
                
                size_t index = 0;
                for (size_t row = (topN == NOBODY ? 0 : 1); row < height - (bottomN == NOBODY ? 0 : 1); ++row) {
                    for (size_t col = 0; col < width; ++col, ++index) {
                        *mfout << curr[index];
                        if (col < width - 1) {
                            *mfout << " ";
                        }
                    }
                    *mfout << "\n";
                }
                if (myCoord == numProcs - 1) {
                    *mfout << "\n";
                }

                mfout->close();
                mfout = NULL;
            }
            MPI_Barrier(comm);
        }
    }
}

#pragma mark - Balancing

void Field::cleanWeights() {
    memset(weights, 0, sizeof(double) * std::max(width, height));
}

void Field::syncWeights() {
    
}

bool Field::balanceNeeded() {
    return false;
}

void Field::balance() {

}


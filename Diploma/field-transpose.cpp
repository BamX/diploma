//
//  Copyright (c) 2016 Nikolay Volosatov. All rights reserved.
//

#include "field-transpose.h"
#include "algo.h"
#include "balancing.h"
#include <cmath>

void FieldTranspose::init() {
    Field::init();

    MPI_Datatype mpi_col_type, mpi_tmp_type;
    MPI_Type_vector((int)height, (int)1, (int)width, MPI_DOUBLE, &mpi_tmp_type);
    MPI_Type_create_resized(mpi_tmp_type, 0, (int)1 * sizeof(double), &mpi_col_type);
    MPI_Type_free(&mpi_tmp_type);

    MPI_Type_contiguous((int)height, mpi_col_type, &mpiVertType);
    MPI_Type_free(&mpi_col_type);
    MPI_Type_commit(&mpiVertType);

    MPI_Type_vector((int)height, (int)height, (int)width, MPI_DOUBLE, &mpi_tmp_type);
    MPI_Type_create_resized(mpi_tmp_type, 0, (int)height * sizeof(double), &mpiHorizType);
    MPI_Type_free(&mpi_tmp_type);
    MPI_Type_commit(&mpiHorizType);

    printf("I'm %d(%d)\twith w:%zu\th:%zu w:%zu\th:%zu.\tTop:%d\tbottom:%d\n",
           myId, ::getpid(), width, height, mySX, mySY, topN, bottomN);
}

FieldTranspose::~FieldTranspose() {
    MPI_Type_free(&mpiHorizType);
    MPI_Type_free(&mpiVertType);
    MPI_Comm_free(&balanceComm);
}

void FieldTranspose::calculateNBS() {
    Field::calculateNBS();

    height = ceil((double)width / numProcs);
    width = height * numProcs;
    hX = algo::ftr().X1() / (width - 1);
    hY = algo::ftr().X2() / (height * numProcs - 1);

    mySY = height * myCoord;
    mySX = 0;

    MPI_Comm_dup(comm, &balanceComm);
}

#pragma mark - Logic

void FieldTranspose::transpose() {
    transpose(transposed ? curr : prev);

    std::swap(hX, hY);
    std::swap(mySX, mySY);
    transposed = transposed == false;
}

size_t FieldTranspose::solveRows() {
    size_t maxIterationsCount = 0;

    for (size_t row = 0; row < height; ++row) {
        fillFactors(row, true);
        double delta = solve(row, true);
        size_t iterationsCount = 1;

        while (delta > epsilon) {
            fillFactors(row, false);
            delta = solve(row, false);
            ++iterationsCount;

            if (iterationsCount > MAX_ITTERATIONS_COUNT) {
                break;
            }
        }
        weights[myCoord * height + row] = iterationsCount;
        maxIterationsCount = std::max(maxIterationsCount, iterationsCount);
    }

    return maxIterationsCount;
}

#pragma mark - MPI

void FieldTranspose::transpose(double *arr) {
    memcpy(buff, arr, width * height * sizeof(double));
    MPI_Alltoall(buff, 1, mpiVertType, arr, 1, mpiHorizType, comm);
}

#pragma mark - Print

double FieldTranspose::view(double x1, double x2) {
    ssize_t x1index = floor(x1 / hX) - mySX;
    ssize_t x2index = floor(x2 / hY) - mySY;

    bool notInMyX1 = x1index < 0 || x1index >= width;
    bool notInMyX2 = x2index < 0 || x2index >= height;
    if (notInMyX1 || notInMyX2) {
        return NOTHING;
    }

    return curr[x2index * width + x1index];
}

void FieldTranspose::printConsole() {
    if (algo::ftr().EnableConsole()) {
        double viewValue = Field::view(algo::ftr().DebugView());

        if (fabs(viewValue - NOTHING) > __DBL_EPSILON__) {
            printf("Field[%d] (itrs: %zu, time: %.5f)\tview: %.7f\n",
                   myId, lastIterrationsCount, t, viewValue);
        }
    }
}

void FieldTranspose::printMatrix() {
    if (algo::ftr().EnableMatrix()) {
        for (int p = 0; p < numProcs; ++p) {
            if (p == myCoord) {
                if (mfout != NULL) {
                    mfout->close();
                }
                mfout = new std::ofstream("matrix.csv", std::ios::app);

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

void FieldTranspose::syncWeights() {
    MPI_Allgather(weights + myCoord * height, (int)height, MPI_DOUBLE,
                  weights, (int)height, MPI_DOUBLE, balanceComm);
}

bool FieldTranspose::balanceNeeded() {
    std::vector<size_t> buckets = balancing::partition(weights, height * numProcs, numProcs);
    if (bfout != NULL) {
        for (size_t i = 0; i < numProcs; ++i) {
            *bfout << buckets[i];
            if (i < numProcs - 1) {
                *bfout << ",";
            }
        }
        *bfout << "\n";
        bfout->flush();
    }
    return false;
}

void FieldTranspose::balance() {

}

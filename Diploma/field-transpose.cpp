//
//  Copyright (c) 2016 Nikolay Volosatov. All rights reserved.
//

#include "field-transpose.h"
#include "algo.h"
#include <cmath>

void FieldTranspose::init() {
    Field::init();

    MPI_Datatype mpi_all_unaligned_t;
    MPI_Type_vector((int)height, (int)height, (int)width, MPI_DOUBLE, &mpi_all_unaligned_t);
    MPI_Type_create_resized(mpi_all_unaligned_t, 0, (int)height * sizeof(double), &mpiAllType);
    MPI_Type_free(&mpi_all_unaligned_t);
    MPI_Type_commit(&mpiAllType);

    printf("I'm %d(%d)\twith w:%zu\th:%zu.\tTop:%d\tbottom:%d\n",
           myId, ::getpid(), width, height, topN, bottomN);
}

FieldTranspose::~FieldTranspose() {
    MPI_Type_free(&mpiAllType);
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
        maxIterationsCount = std::max(maxIterationsCount, iterationsCount);
    }

    return maxIterationsCount;
}

#pragma mark - MPI

void FieldTranspose::transpose(double *arr) {
    for (size_t k = 0; k < numProcs; ++k) {
        for (size_t i = 0; i < height; ++i) {
            for (size_t j = 0; j < height; ++j) {
                size_t index = i * width + k * height + j;
                size_t newIndex = j * width + k * height + i;
                buff[newIndex] = arr[index];
            }
        }
    }
    MPI_Alltoall(buff, 1, mpiAllType, arr, 1, mpiAllType, comm);
}

#pragma mark - Print

void FieldTranspose::printConsole() {
    if (algo::ftr().EnableConsole()) {
        double viewValue = view(algo::ftr().DebugView());

        if (fabs(viewValue - NOTHING) > __DBL_EPSILON__) {
            printf("Field[%d] (itrs: %zu, time: %.5f)\tview: %.7f\n",
                   myId, lastIterrationsCount, t, viewValue);
        }
    }
}

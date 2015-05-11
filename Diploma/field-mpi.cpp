//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#include "field.h"
#include <cmath>
#include <sys/types.h>
#include <unistd.h>

//#define DEBUG_PRINT
//#define DEBUG_WAIT


void Field::debug(const char *name) {
#ifdef DEBUG_PRINT
    printf("I'm %d before %s\n", myId, name);
#endif
}

void Field::calculateNBS() {
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    int dims[] = { numProcs };
    int wrap[] = { 0 };
    MPI_Cart_create(MPI_COMM_WORLD, 1, dims, wrap, 1, &comm);

    MPI_Comm_rank(comm, &myId);
    MPI_Cart_coords(comm, myId, 1, &myCoord);
    MPI_Cart_shift(comm, 0, 1, &topN, &bottomN);

    // Workaround for equal new heights
    height = ceil((double)height / numProcs) * numProcs;
    width = height;
    hX = ftr.X1() / (width - 1);
    hY = ftr.X2() / (height - 1);

    size_t newHeight = height / numProcs;

    mySY = newHeight * myCoord;
    mySX = 0;

    if (bottomN == NOBODY) {
        newHeight = height - newHeight * (numProcs - 1);
    }

    height = newHeight;

    MPI_Datatype mpi_all_unaligned_t;
    MPI_Type_vector((int)height, (int)height, (int)width, MPI_DOUBLE, &mpi_all_unaligned_t);
    MPI_Type_create_resized(mpi_all_unaligned_t, 0, (int)height * sizeof(double), &mpiAllType);
    MPI_Type_free(&mpi_all_unaligned_t);
    MPI_Type_commit(&mpiAllType);

    printf("I'm %d(%d)\twith w:%zu\th:%zu.\tTop:%d\tbottom:%d\n",
           myId, ::getpid(), width, height, topN, bottomN);

#ifdef DEBUG_WAIT
    int waiter = myId;
    while (waiter == WAITER) sleep(5);
    printf("GO\n");
#endif
}

void Field::finalize() {
    MPI_Type_free(&mpiAllType);
}

void Field::transpose(double *arr) {
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

void Field::reduceViews() {
    for (size_t index = 0, len = ftr.ViewCount(); index < len; ++index) {
        double value = view(index);
        double result = 0;
        MPI_Reduce(&value, &result, 1, MPI_DOUBLE, MPI_MAX, MASTER, comm);
        views[index] = result;
    }
}

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

void Field::initFactors() {
    int procs, id;
    MPI_Comm_size(MPI_COMM_WORLD, &procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    for (size_t p = 0; p < procs; ++p) {
        if (p == id) {
            ftr.initFactors();
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

void Field::printMatrix() {
    if (ftr.EnableMatrix()) {
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

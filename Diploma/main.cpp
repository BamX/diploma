//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#include <stdio.h>
#include <mpi.h>

#include "field-static.h"
#include "field-transpose.h"
#include "factors.h"

void createVType(size_t width, size_t height, size_t bWidth, MPI_Datatype *type) {
    MPI_Datatype mpi_retmp_type;
    MPI_Datatype mpi_col_type;

    MPI_Type_vector((int)height, (int)1, (int)width, MPI_DOUBLE, &mpi_retmp_type);

    MPI_Type_create_resized(mpi_retmp_type, 0, (int)1 * sizeof(double), &mpi_col_type);
    MPI_Type_free(&mpi_retmp_type);

    MPI_Type_contiguous((int)bWidth, mpi_col_type, type);
    MPI_Type_commit(type);
    MPI_Type_free(&mpi_col_type);
}

void createHType(size_t width, size_t height, size_t bWidth, MPI_Datatype *type) {
    MPI_Datatype mpi_tmp_type;

    MPI_Type_vector((int)height, (int)bWidth, (int)width, MPI_DOUBLE, &mpi_tmp_type);

    MPI_Type_create_resized(mpi_tmp_type, 0, (int)height * sizeof(double), type);
    MPI_Type_commit(type);
    MPI_Type_free(&mpi_tmp_type);
}

int __main(int argc, char * argv[]) {
    MPI_Init(&argc, &argv);

    const int numProc = 4;
    const int width = 12;
    const int height = width / numProc;

    int myId = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myId);

    double arr[width * width * 100] = {0}, barr[width * width * 100] = {0};
    for (int i = 0; i < height * width; ++i) {
        arr[i] = height * width * myId + i;
    }

    int sendcounts[numProc] = {1};
    int recvcounts[numProc] = {1};
    int senddispls[numProc] = {1};
    int recvdispls[numProc] = {1};
    MPI_Datatype sendtypes[numProc], recvtypes[numProc];

    int hBuckets[numProc] = { 3, 2, 4, 3 };
    int vBuckets[numProc] = { 3, 3, 3, 3 };

    for (size_t i = 0; i < numProc; ++i) {
        sendcounts[i] = recvcounts[i] = 1;
        senddispls[i] = i == 0 ? 0 : (int)(senddispls[i - 1] + hBuckets[i - 1] * sizeof(double));
        recvdispls[i] = i == 0 ? 0 : (int)(recvdispls[i - 1] + vBuckets[i - 1] * sizeof(double));
        createVType(width, vBuckets[myId], hBuckets[i], sendtypes + i);
        createHType(width, hBuckets[myId], vBuckets[i], recvtypes + i);
        std::cerr << "PROC " << myId << " V:" << vBuckets[myId] << "x" << hBuckets[i] << " H:" << vBuckets[myId] << "x" << hBuckets[i] << "\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);

    for (int k = 0; k < numProc; ++k) {
        if (k == myId) {
            for (int i = 0; i < height; ++i) {
                for (int j = 0; j < width; ++j) {
                    std::cout << arr[i * width + j] << "\t";
                }
                std::cout << "\n";
            }
            std::cout << "\n";
        }
        std::cout.flush();
        MPI_Barrier(MPI_COMM_WORLD);
    }
    if (myId == 0) {
        std::cout << "\n\n";
    }

    memcpy(barr, arr, width * height * sizeof(double));
    MPI_Alltoallw(barr, sendcounts, senddispls, sendtypes, arr, recvcounts, recvdispls, recvtypes, MPI_COMM_WORLD);

    for (int k = 0; k < numProc; ++k) {
        if (k == myId) {
            for (int i = 0; i < hBuckets[k]; ++i) {
                for (int j = 0; j < width; ++j) {
                    std::cout << arr[i * width + j] << "\t";
                }
                std::cout << "\n";
            }
            std::cout << "\n";
        }
        std::cout.flush();
        MPI_Barrier(MPI_COMM_WORLD);
    }
    if (myId == 0) {
        std::cout << "\n\n";
    }

    MPI_Finalize();
    return 0;
}

int main(int argc, char * argv[]) {
    MPI_Init(&argc, &argv);

    {
        FieldTranspose field;
        //FieldStatic field;
        
        field.init();

        {
            const auto startTime = std::clock();

            while (field.done() == false) {
                field.solve();
            }

            const auto endTime = std::clock();

            int rank = 0;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            if (rank == 0) {
                std::cerr << "time: " << double(endTime - startTime) / CLOCKS_PER_SEC << '\n';
            }
        }

        field.finalize();
    }
    
    MPI_Finalize();
    return 0;
}

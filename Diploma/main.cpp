//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#include <stdio.h>
#include <mpi.h>

#include "field-static.h"
#include "field-transpose.h"
#include "factors.h"
#include "algo.h"

#ifndef ALGV
#define ALGV 0
#endif

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

void initConfig(int argc, char * argv[]) {
    char configFilename[100] = {0};
    if (argc >= 2) {
        strcpy(configFilename, argv[1]);
    } else {
        strcpy(configFilename, "config.ini");
    }

    int procs, id;
    MPI_Comm_size(MPI_COMM_WORLD, &procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    for (size_t p = 0; p < procs; ++p) {
        if (p == id) {
            algo::ftr().initFactors(configFilename);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

int main(int argc, char * argv[]) {
    MPI_Init(&argc, &argv);

    initConfig(argc, argv);

    auto startTime = bx_clock_t::now();
    {
        Field *field = NULL;
        if (algo::ftr().Algorithm() == kAlgorithmTranspose) {
            field = new FieldTranspose;
        } else if (algo::ftr().Algorithm() == kAlgorithmStatic) {
            field = new FieldStatic;
        }

        for (size_t k = 0; k < algo::ftr().Repeats(); ++k) {
            field->init();

            while (field->done() == false) {
                field->solve();
            }

            field->finalize();
        }

        delete field;
    }

    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Barrier(MPI_COMM_WORLD);
    auto picosecCount = std::chrono::duration<unsigned long long, std::pico>(bx_clock_t::now() - startTime).count();
    if (myRank == 0) {
        std::cerr << picosecCount * 1.e-12 << "\n";
    }

    MPI_Finalize();
    return 0;
}

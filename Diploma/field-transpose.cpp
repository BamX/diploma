//
//  Copyright (c) 2016 Nikolay Volosatov. All rights reserved.
//

#include "field-transpose.h"
#include "algo.h"
#include "balancing.h"
#include <cmath>

#include <sys/types.h>
#include <unistd.h>

static int const kDefaultBalancingCounter = 16;

void FieldTranspose::init() {
    Field::init();

    balancingCounter = kDefaultBalancingCounter;

    printf("I'm %d(%d)\twith w:%zu\th:%zu w:%zu\th:%zu.\tTop:%d\tbottom:%d\n",
           myId, ::getpid(), width, height, mySX, mySY, topN, bottomN);
}

FieldTranspose::~FieldTranspose() {
    delete[] sendcounts;
    delete[] recvcounts;
    delete[] senddispls;
    delete[] recvdispls;
    delete[] gathercounts;
    delete[] gatherdispls;
    for (size_t i = 0; i < numProcs; ++i) {
        MPI_Type_free(sendtypes + i);
        MPI_Type_free(recvtypes + i);
    }
    delete[] sendtypes;
    delete[] recvtypes;

    delete[] weights;
    delete[] weightsT;

    MPI_Comm_free(&balanceComm);
}

void FieldTranspose::calculateNBS() {
    Field::calculateNBS();

    height = ceil((double)width / numProcs);
    width = height * numProcs;
    hX = algo::ftr().X1() / (width - 1);
    hY = algo::ftr().X2() / (height * numProcs - 1);

    mySY = mySYT = height * myCoord;
    mySX = 0;

    heightCapacity = height;

    hBuckets = new size_t[numProcs];
    vBuckets = new size_t[numProcs];
    nextBuckets.resize(numProcs);
    nextBucketsT.resize(numProcs);
    for (size_t i = 0; i < numProcs; ++i) {
        hBuckets[i] = vBuckets[i] = nextBuckets[i] = nextBucketsT[i] = (int)height;
    }

    sendcounts = new int[numProcs];
    senddispls = new int[numProcs];
    recvcounts = new int[numProcs];
    recvdispls = new int[numProcs];
    gathercounts = new int[numProcs];
    gatherdispls = new int[numProcs];
    sendtypes = new MPI_Datatype[numProcs];
    recvtypes = new MPI_Datatype[numProcs];

    for (size_t i = 0; i < numProcs; ++i) {
        sendcounts[i] = recvcounts[i] = 1;
        senddispls[i] = i == 0 ? 0 : (int)(senddispls[i - 1] + hBuckets[i - 1] * sizeof(double));
        recvdispls[i] = i == 0 ? 0 : (int)(recvdispls[i - 1] + vBuckets[i - 1] * sizeof(double));

        createVType(width, vBuckets[myCoord], hBuckets[i], sendtypes + i);
        createHType(width, hBuckets[myCoord], vBuckets[i], recvtypes + i);
    }

    weights = new double[width];
    memset(weights, 0, width * sizeof(double));
    weightsT = new double[width];
    memset(weightsT, 0, width * sizeof(double));

    MPI_Comm_dup(comm, &balanceComm);
}

#pragma mark - Logic

void FieldTranspose::transpose() {
    transpose(transposed ? curr : prev);

    std::swap(hX, hY);
    std::swap(mySY, mySYT);
    std::swap(hBuckets, vBuckets);
    std::swap(senddispls, recvdispls);
    std::swap(sendtypes, recvtypes);
    std::swap(nextBuckets, nextBucketsT);
    std::swap(weights, weightsT);

    transposed = transposed == false;
}

size_t FieldTranspose::solveRows() {
    debug(0).flush();
    
    size_t maxIterationsCount = 0;

    for (size_t row = 0; row < height; ++row) {
        fillFactors(row, true);
        double delta = solve(row, true);
        size_t iterationsCount = 1;
        auto startTime = picosecFromStart();

        while (delta > epsilon) {
            fillFactors(row, false);
            delta = solve(row, false);
            ++iterationsCount;

            if (iterationsCount > MAX_ITTERATIONS_COUNT) {
                break;
            }
        }
        //debug() << "Write to " << mySY + row << " of " << width << "\n";
        weights[mySY + row] += iterationsCount;// + (picosecFromStart() - startTime) * 1e-12 / iterationsCount;
        maxIterationsCount = std::max(maxIterationsCount, iterationsCount);
    }

    return maxIterationsCount;
}

#pragma mark - Print

double FieldTranspose::view(double x1, double x2) {
    long x1index = floor(x1 / hX) - mySX;
    long x2index = floor(x2 / hY) - mySY;

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
            printf("Field[%d] (itrs: %zu, time: %.5f) ctime: %.1f\tview: %.7f\n",
                   myId, lastIterrationsCount, t, picosecFromStart() * 1e-12, viewValue);
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
                    *mfout << " @" << myCoord << "\n";
                }
                if (myCoord == numProcs - 1) {
                    *mfout << "\n";
                }

                mfout->close();
                delete mfout;
                mfout = NULL;
            }
            MPI_Barrier(comm);
        }
    }
}

#pragma mark - MPI

void FieldTranspose::transpose(double *arr) {
    //resize(hBuckets[myCoord]);
    memcpy(buff, arr, height * width * sizeof(double));
    height = hBuckets[myCoord];

    //debug() << "Alltoallw\n";

    /*for (size_t i = 0; i < numProcs; ++i) {
        int sendsize, recvsize;
        char sendbuff[255], recvbuff[255];
        MPI_Type_get_name(sendtypes[i], sendbuff, &sendsize);
        MPI_Type_get_name(recvtypes[i], recvbuff, &recvsize);
        MPI_Type_size(sendtypes[i], &sendsize);
        MPI_Type_size(recvtypes[i], &recvsize);

        debug() << "to: " << i << " count: " << sendcounts[i] << " disp: " << senddispls[i] << " size: " << sendsize / 8 << " name: " << sendbuff
        <<" | " << "from: " << i << " count: " << recvcounts[i] << " disp: " << recvdispls[i] << " size: " << recvsize / 8  << " name: " << recvbuff << "\n";
    }*/

    MPI_Alltoallw(buff, sendcounts, senddispls, sendtypes, arr, recvcounts, recvdispls, recvtypes, comm);

    //debug() << "OK height: " << height << "\n";
}

#pragma mark - Balancing

void FieldTranspose::resize(size_t newHeight) {
    if (newHeight <= heightCapacity) {
        return;
    }

    while (heightCapacity < newHeight) {
        heightCapacity = heightCapacity * 10;
    }

    memcpy(buff, prev, width * height * sizeof(double));
    delete[] prev;
    prev = new double[heightCapacity * width];
    memcpy(prev, buff, width * height * sizeof(double));

    memcpy(buff, curr, width * height * sizeof(double));
    delete[] curr;
    curr = new double[heightCapacity * width];
    memcpy(curr, buff, width * height * sizeof(double));

    delete[] buff;
    buff = new double[heightCapacity * width];

    delete[] maF;
    delete[] mbF;
    delete[] mcF;
    delete[] mfF;

    maF = new double[heightCapacity * width];
    mbF = new double[heightCapacity * width];
    mcF = new double[heightCapacity * width];
    mfF = new double[heightCapacity * width];

    //debug() << "NEW HEIGHT " << heightCapacity << "\n";
}

void FieldTranspose::createVType(size_t width, size_t height, size_t bWidth, MPI_Datatype *type) {
    MPI_Datatype mpi_retmp_type;
    MPI_Datatype mpi_col_type;

    MPI_Type_vector((int)height, (int)1, (int)width, MPI_DOUBLE, &mpi_retmp_type);

    MPI_Type_create_resized(mpi_retmp_type, 0, (int)1 * sizeof(double), &mpi_col_type);
    MPI_Type_free(&mpi_retmp_type);

    MPI_Type_contiguous((int)bWidth, mpi_col_type, type);
    //MPI_Type_set_name(*type, "V_Block");
    MPI_Type_commit(type);
    MPI_Type_free(&mpi_col_type);
}

void FieldTranspose::createHType(size_t width, size_t height, size_t bWidth, MPI_Datatype *type) {
    MPI_Datatype mpi_tmp_type;

    MPI_Type_vector((int)height, (int)bWidth, (int)width, MPI_DOUBLE, &mpi_tmp_type);

    MPI_Type_create_resized(mpi_tmp_type, 0, (int)height * sizeof(double), type);
    //MPI_Type_set_name(*type, "H_Block");
    MPI_Type_commit(type);
    MPI_Type_free(&mpi_tmp_type);
}

void FieldTranspose::syncWeights() {
    if (algo::ftr().Balancing()) {
        //debug() << "SW ..." << " " << myId << " " << mySY << " " << height << "\n";
        for (size_t i = 0; i < numProcs; ++i) {
            gathercounts[i] = (int)vBuckets[i];
            gatherdispls[i] = i == 0 ? 0 : (int)(gatherdispls[i - 1] + gathercounts[i - 1]);
            //debug() << "SWW ..." << " " << myId << " " << i << " " << height << " " << gathercounts[i] << " " << gatherdispls[i] << "\n";
        }
        
        MPI_Gatherv(weights + mySY, (int)height, MPI_DOUBLE, weights, gathercounts, gatherdispls, MPI_DOUBLE, MASTER, balanceComm);

        if (myId == MASTER) {
            nextBucketsT = balancing::partition(weights, width, numProcs);
        }

        MPI_Bcast(&nextBucketsT[0], (int)numProcs, MPI_INT, MASTER, balanceComm);

        //debug() << "weights: ";
        for (int i = 0; i < width; ++i) {
            debug(0) << weights[i] << " ";
        }
        //debug(0) << "\n";

        //debug() << "buckets: ";
        for (int i = 0; i < numProcs; ++i) {
            debug(0) << nextBucketsT[i] << " ";
        }
        //debug(0) << "\n";
        //debug() << "SW OK\n";

        memset(weights, 0, std::max(height, width) * sizeof(double));
    } else {
        for (size_t i = 0; i < numProcs; ++i) {
            nextBucketsT[i] = (int)(width / numProcs);
        }
    }
}

bool FieldTranspose::balanceNeeded() {
    if (algo::ftr().Balancing()) {
        if (bfout != NULL) {
            for (size_t i = 0; i < numProcs; ++i) {
                *bfout << nextBuckets[i];
                if (i < numProcs - 1) {
                    *bfout << ",";
                }
            }
            *bfout << "\n";
            bfout->flush();
        }

        balancingCounter -= 1;
        if (balancingCounter < 0) {
           balancingCounter = kDefaultBalancingCounter;
        }
        return balancingCounter <= 1;
    } else {
        return false;
    }
}

void FieldTranspose::balance() {
    if (myCoord == 0) {
        //debug() << "=================\n";
    }

    //debug() << "set-buckets: ";
    for (int i = 0; i < numProcs; ++i) {
        //debug(0) << nextBuckets[i] << " ";
    }
    //debug(0) << "\n";
    //debug() << "SW OK\n";

    hBuckets[myCoord] = nextBuckets[myCoord];
    mySYT = 0;
    for (size_t i = 0; i < numProcs; ++i) {
        hBuckets[i] = nextBuckets[i];

        if (i < myCoord) {
            mySYT += hBuckets[i];
        }
        
        senddispls[i] = i == 0 ? 0 : (int)(senddispls[i - 1] + hBuckets[i - 1] * sizeof(double));
        recvdispls[i] = i == 0 ? 0 : (int)(recvdispls[i - 1] + vBuckets[i - 1] * sizeof(double));

        MPI_Type_free(sendtypes + i);
        MPI_Type_free(recvtypes + i);

        createVType(width, vBuckets[myCoord], hBuckets[i], sendtypes + i);
        createHType(width, hBuckets[myCoord], vBuckets[i], recvtypes + i);

        /*debug() << "PROC " << myCoord << " V:" << vBuckets[myCoord] << "x" << hBuckets[i]
                                      << " H:" << hBuckets[myCoord] << "x" << vBuckets[i]
                                      << "  " << mySYT << "\n";
         */
    }
}

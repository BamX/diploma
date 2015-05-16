//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#include "field.h"
#include <cmath>
#include <sys/types.h>
#include <unistd.h>

//#define DEBUG_PRINT
//#define DEBUG_WAIT

static int const TAG_PREFIX = 1000000;
static int const TAG_TOP_TO_BOTTOM_ROWS = 21;
static int const TAG_BOTTOM_TO_TOP_ROWS = 22;
static int const TAG_CALCULATING_ROWS = 23;
static int const TAG_FIRST_PASS = 24 * TAG_PREFIX;
static int const TAG_SECOND_PASS = 25 * TAG_PREFIX;

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
    rightN = leftN = NOBODY;

    size_t newHeight = height / numProcs;

    mySY = newHeight * myCoord;
    mySX = 0;

    if (bottomN == NOBODY) {
        newHeight = height - newHeight * (numProcs - 1);
    }
    newHeight += (topN != NOBODY ? 1 : 0) + (bottomN != NOBODY ? 1 : 0);

    height = newHeight;
    bundleSizeLimit = std::max(ceil((double)width / numProcs / 2), 15.0);

    printf("I'm %d(%d)\twith w:%zu\th:%zu\tbs:%zu.\tTop:%d\tbottom:%d\n",
           myId, ::getpid(), width, height, bundleSizeLimit, topN, bottomN);

#ifdef DEBUG_WAIT
    int waiter = myId;
    while (waiter == WAITER) sleep(5);
    printf("GO\n");
#endif
}

void Field::finalize() {
}

void Field::sendRecieveCalculatingRows() {
    MPI_Recv(calculatingRows, (int)height, MPI::BOOL, leftN, TAG_CALCULATING_ROWS, comm, MPI_STATUS_IGNORE);
    MPI_Request request;
    MPI_Isend(calculatingRows, (int)height, MPI::BOOL, rightN, TAG_CALCULATING_ROWS, comm, &request);
}

void Field::sendFistPass(size_t fromRow) {
    // calculatingRows x [prevCalculatingRows] + (y + b + c + f) x [calculatingRows]
    if (rightN != NOBODY) {
        double *sBuff = sendBuff + fromRow * SEND_PACK_SIZE;

        int sSize = 0;
        size_t bundleSize = 0;
        for (size_t row = fromRow; row < height; ++row) {
            if (calculatingRows[row] == false) {
                continue;
            }

            size_t index = row * width + width - 2;
            sBuff[sSize++] = mbF[index];
            sBuff[sSize++] = mcF[index];
            sBuff[sSize++] = mfF[index];
            ++bundleSize;
            if (bundleSize >= bundleSizeLimit) {
                break;
            }
        }

        if (sSize == 0) {
            return;
        }
        MPI_Request request;
        MPI_Isend(sBuff, sSize, MPI_DOUBLE, rightN, TAG_FIRST_PASS + (int)fromRow, comm, &request);
    }
}

void Field::recieveFirstPass(size_t fromRow) {
    // calculatingRows x [prevCalculatingRows] + (y + b + c + f) x [calculatingRows]
    if (leftN != NOBODY) {
        int sSize = 0;
        size_t bundleSize = 0;
        for (size_t row = fromRow; row < height; ++row) {
            if (calculatingRows[row] == false) {
                continue;
            }
            sSize += 3;
            ++bundleSize;

            if (bundleSize >= bundleSizeLimit) {
                break;
            }
        }
        if (sSize == 0) {
            return;
        }

        MPI_Recv(receiveBuff, sSize, MPI_DOUBLE, leftN, TAG_FIRST_PASS + (int)fromRow, comm, MPI_STATUS_IGNORE);
        sSize = 0;
        bundleSize = 0;
        for (size_t row = fromRow; row < height; ++row) {
            if (calculatingRows[row] == false) {
                continue;
            }

            size_t index = row * width;
            mbF[index] = receiveBuff[sSize++];
            mcF[index] = receiveBuff[sSize++];
            mfF[index] = receiveBuff[sSize++];
            ++bundleSize;
            if (bundleSize >= bundleSizeLimit) {
                break;
            }
        }
    }
}

void Field::sendSecondPass(size_t fromRow) {
    // (calculatingRows + y) x [prevCalculatingRows]
    if (leftN != NOBODY) {
        double *sBuff = sendBuff + fromRow * SEND_PACK_SIZE;
        int sSize = 0;
        size_t bundleSize = 0;
        for (size_t row = fromRow; row < height; ++row) {
            if (prevCalculatingRows[row] == false) {
                continue;
            }
            sBuff[sSize++] = (calculatingRows[row] ? 1 : -1) * curr[row * width + 1];
            ++bundleSize;
            if (bundleSize >= bundleSizeLimit) {
                break;
            }
        }
        if (sSize == 0) {
            return;
        }

        MPI_Request request;
        MPI_Isend(sBuff, sSize, MPI_DOUBLE, leftN, TAG_SECOND_PASS + (int)fromRow, comm, &request);
    }
}

void Field::recieveSecondPass(size_t fromRow) {
    // (calculatingRows + y) x [prevCalculatingRows]
    if (rightN != NOBODY) {
        int sSize = 0;
        size_t bundleSize = 0;
        for (size_t row = fromRow; row < height; ++row) {
            if (prevCalculatingRows[row] == false) {
                continue;
            }
            ++sSize;
            ++bundleSize;
            if (bundleSize >= bundleSizeLimit) {
                break;
            }
        }
        if (sSize == 0) {
            return;
        }

        MPI_Recv(receiveBuff, sSize, MPI_DOUBLE, rightN, TAG_SECOND_PASS + (int)fromRow, comm, MPI_STATUS_IGNORE);
        sSize = 0;
        bundleSize = 0;
        for (size_t row = fromRow; row < height; ++row) {
            if (prevCalculatingRows[row] == false) {
                continue;
            }

            double value = receiveBuff[sSize++];
            calculatingRows[row] = value > 0;
            curr[(row + 1) * width - 1] = calculatingRows[row] ? value : -value;
            ++bundleSize;
            if (bundleSize >= bundleSizeLimit) {
                break;
            }
        }
    }
}

void Field::sendRecieveRows() {
    MPI_Sendrecv(prev + (height - 2) * width, (int)width, MPI_DOUBLE, bottomN, TAG_TOP_TO_BOTTOM_ROWS,
                 prev, (int)width, MPI_DOUBLE, topN, TAG_TOP_TO_BOTTOM_ROWS,
                 comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(prev + width, (int)width, MPI_DOUBLE, topN, TAG_BOTTOM_TO_TOP_ROWS,
                 prev + (height - 1) * width, (int)width, MPI_DOUBLE, bottomN, TAG_BOTTOM_TO_TOP_ROWS,
                 comm, MPI_STATUS_IGNORE);
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
        if (t > nextFrameTime) {
            nextFrameTime += ftr.totalTime() / ftr.MatrixFramesCount();

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
}

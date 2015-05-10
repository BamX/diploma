//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#include "field.h"
#include <cmath>
#include <sys/types.h>
#include <unistd.h>

//#define DEBUG_PRINT
//#define DEBUG_WAIT

static int const MAX_ROW_INDEX = 1000000;
static int const TAG_TOP_TO_BOTTOM_ROWS = 23;
static int const TAG_BOTTOM_TO_TOP_ROWS = 24;
static int const TAG_LEFT_TO_RIGHT_ROW = 1 * MAX_ROW_INDEX;
static int const TAG_FIRST_PASS = 3 * MAX_ROW_INDEX;
static int const TAG_SECOND_PASS = 4 * MAX_ROW_INDEX;

void deltasReducer(void *in, void *inout, int *len, MPI_Datatype *dptr) {
    double *din = (double *)in;
    double *dinout = (double *)inout;

    for (int i = 0; i < *len; ++i) {
        *dinout = std::max(*din, *dinout);
        ++din; ++dinout;
    }
}

void Field::debug(const char *name) {
#ifdef DEBUG_PRINT
    printf("I'm %d before %s\n", myId, name);
#endif
}

void Field::calculateNBS() {
    MPI_Op_create(&deltasReducer, 1, &deltasReducerOp);

    int numProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    int stX = 0, stY = 0;
    calculateGrid(numProcs, (double)width / height, stX, stY);

    int dims[2] = { stY, stX };
    int wrap[2] = { 0, 0 };
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, wrap, 1, &comm);

    int coord[2];
    MPI_Comm_rank(MPI_COMM_WORLD, &myId);
    MPI_Cart_coords(comm, myId, 2, coord);
    MPI_Cart_shift(comm, 0, 1, &topN, &bottomN);
    MPI_Cart_shift(comm, 1, 1, &leftN, &rightN);

    createRowColComms(coord[0], coord[1], stX, stY);

    size_t newHeight = height / stY;
    size_t newWidth = width / stX;

    mySY = newHeight * coord[0];
    mySX = newWidth * coord[1];

    if (bottomN == NOBODY) {
        newHeight = height - newHeight * (stY - 1);
    }
    newHeight += (topN != NOBODY ? 1 : 0) + (bottomN != NOBODY ? 1 : 0);

    if (rightN == NOBODY) {
        newWidth = width - newWidth * (stX - 1);
    }
    newWidth += (rightN != NOBODY ? 1 : 0) + (leftN != NOBODY ? 1 : 0);

    width = newWidth;
    height = newHeight;

    printf("I'm %d(%d)\twith w:%zu\th:%zu.\tTop:%d\tbottom:%d\tleft:%d\tright:%d\n",
           myId, ::getpid(), width, height, topN, bottomN, leftN, rightN);

#ifdef DEBUG_WAIT
    int waiter = myId;
    while (waiter == WAITER) sleep(5);
    MPI_Barrier(comm);
    printf("GO\n");
#endif
}

void Field::createRowColComms(int myI, int myJ, int stX, int stY) {
    MPI_Group commGroup;
    MPI_Comm_group(comm, &commGroup);

    int *colRanks = new int[stY];
    for (int i = 0; i < stY; ++i) {
        int coord[] = { i, myJ };
        MPI_Cart_rank(comm, coord, colRanks + i);
    }
    MPI_Group colGroup;
    MPI_Group_incl(commGroup, stY, colRanks, &colGroup);
    MPI_Comm_create(comm, colGroup, &colComm);

    int *rowRanks = new int[stX];
    for (int j = 0; j < stX; ++j) {
        int coord[] = { myI, j };
        MPI_Cart_rank(comm, coord, rowRanks + j);
    }
    MPI_Group rowGroup;
    MPI_Group_incl(commGroup, stX, rowRanks, &rowGroup);
    MPI_Comm_create(comm, rowGroup, &rowComm);
}

void Field::calculateGrid(int numProcs, double stExpected, int &stX, int &stY) {
    stX = 1; stY = numProcs;

    bool withSwap = false;
    if (stExpected > 1) {
        stExpected = 1.0 / stExpected;
        withSwap = true;
    }

    double stFactor = 1.0 / numProcs;
    for (int k = 2, len = numProcs / 2; k <= len; ++k) {
        int y = numProcs / k;
        double newStFactor = (double)k / y;
        if (numProcs % k == 0 && fabs(newStFactor - stExpected) < fabs(stFactor - stExpected)) {
            stX = k; stY = y;
            stFactor = newStFactor;
        }
    }
    
    if (withSwap) {
        std::swap(stX, stY);
    }
}

void Field::sendReceivePrevRows() {
    debug("prevs");
    MPI_Sendrecv(prev + (height - 2) * width, (int)width, MPI_DOUBLE, bottomN, TAG_TOP_TO_BOTTOM_ROWS,
                 prev, (int)width, MPI_DOUBLE, topN, TAG_TOP_TO_BOTTOM_ROWS,
                 comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(prev + width, (int)width, MPI_DOUBLE, topN, TAG_BOTTOM_TO_TOP_ROWS,
                 prev + (height - 1) * width, (int)width, MPI_DOUBLE, bottomN, TAG_BOTTOM_TO_TOP_ROWS,
                 comm, MPI_STATUS_IGNORE);
    debug("prevs");
}

void Field::sendReceiveLeftBorders() {
    size_t shift = 0;
    for (size_t row = (topN != NOBODY ? 1 : 0), len = height - (bottomN != NOBODY ? 1 : 0); row < len; ++row) {
        double *rw = curr + row * width;

        sendBuff[shift++] = rw[width - 2];
    }

    MPI_Sendrecv(sendBuff, (int)shift, MPI_DOUBLE, rightN, TAG_LEFT_TO_RIGHT_ROW,
                 recvBuff, (int)shift, MPI_DOUBLE, leftN, TAG_LEFT_TO_RIGHT_ROW,
                 comm, MPI_STATUS_IGNORE);

    shift = 0;
    for (size_t row = (topN != NOBODY ? 1 : 0), len = height - (bottomN != NOBODY ? 1 : 0); row < len; ++row) {
        double *rw = curr + row * width;

        rw[0] = recvBuff[shift++];
    }
}

void Field::sendFirstPass() {
    if (rightN != NOBODY) {
        size_t shift = 0;
        for (size_t row = (topN != NOBODY ? 1 : 0), len = height - (bottomN != NOBODY ? 1 : 0); row < len; ++row) {

            double *rbF = bF + row * width;
            double *rcF = cF + row * width;
            double *rfF = fF + row * width;

            sendBuff[shift++] = rbF[width - 2];
            sendBuff[shift++] = rcF[width - 2];
            sendBuff[shift++] = rfF[width - 2];
        }

        debug("first pass send");
        MPI_Send(sendBuff, (int)shift, MPI_DOUBLE, rightN, TAG_FIRST_PASS, comm);
        debug("first pass send");
    }
}

void Field::receiveFirstPass() {
    if (leftN != NOBODY) {
        size_t shift = 0;
        for (size_t row = (topN != NOBODY ? 1 : 0), len = height - (bottomN != NOBODY ? 1 : 0); row < len; ++row) {
            shift += 3;
        }

        debug("first pass receive");
        MPI_Recv(recvBuff, (int)shift, MPI_DOUBLE, leftN, TAG_FIRST_PASS, comm, MPI_STATUS_IGNORE);
        debug("first pass receive");

        shift = 0;
        for (size_t row = (topN != NOBODY ? 1 : 0), len = height - (bottomN != NOBODY ? 1 : 0); row < len; ++row) {

            double *rbF = bF + row * width;
            double *rcF = cF + row * width;
            double *rfF = fF + row * width;

            rbF[0] = recvBuff[shift++];
            rcF[0] = recvBuff[shift++];
            rfF[0] = recvBuff[shift++];
        }
    }
}

void Field::sendSecondPass() {
    if (leftN != NOBODY) {
        size_t shift = 0;
        for (size_t row = (topN != NOBODY ? 1 : 0), len = height - (bottomN != NOBODY ? 1 : 0); row < len; ++row) {

            double *y = curr + row * width;

            sendBuff[shift++] = y[1];
        }

        debug("second pass send");
        MPI_Send(sendBuff, (int)shift, MPI_DOUBLE, leftN, TAG_SECOND_PASS, comm);
        debug("second pass send");
    }
}

void Field::receiveSecondPass() {
    if (rightN != NOBODY) {
        size_t shift = 0;
        for (size_t row = (topN != NOBODY ? 1 : 0), len = height - (bottomN != NOBODY ? 1 : 0); row < len; ++row) {
            ++shift;
        }

        debug("second pass receive");
        MPI_Recv(recvBuff, (int)shift, MPI_DOUBLE, rightN, TAG_SECOND_PASS, comm, MPI_STATUS_IGNORE);
        debug("second pass receive");

        shift = 0;
        for (size_t row = (topN != NOBODY ? 1 : 0), len = height - (bottomN != NOBODY ? 1 : 0); row < len; ++row) {

            double *y = curr + row * width;

            y[width - 1] = recvBuff[shift++];
        }
    }
}

double Field::reduceMaxDelta() {
    int count = (int)std::max(width, height);
    memcpy(sendBuff, rowDeltas, count * sizeof(double));
    MPI_Allreduce(sendBuff, rowDeltas, count, MPI_DOUBLE, deltasReducerOp, rowComm);

    double maxDelta = 0;
    for (size_t row = (topN != NOBODY ? 1 : 0), len = height - (bottomN != NOBODY ? 1 : 0); row < len; ++row) {
        maxDelta = std::max(rowDeltas[row], maxDelta);
    }

    //printf("max delta: %.5f\n", maxDelta);
    return maxDelta;
}

void Field::reduceViews() {
    for (size_t index = 0, len = ftr.ViewCount(); index < len; ++index) {
        double value = view(index);
        double result = 0;
        MPI_Reduce(&value, &result, 1, MPI_DOUBLE, MPI_MAX, MASTER, comm);
        views[index] = result;
    }
}

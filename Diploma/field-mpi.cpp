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
static int const TAG_RIGHT_TO_LEFT_ROW = 2 * MAX_ROW_INDEX;
static int const TAG_FIRST_PASS = 3 * MAX_ROW_INDEX;
static int const TAG_SECOND_PASS = 4 * MAX_ROW_INDEX;

void Field::debug(const char *name) {
#ifdef DEBUG_PRINT
    printf("I'm %d before %s\n", myId, name);
#endif
}

void Field::calculateNBS() {
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

void Field::sendReceiveCurrRowLeftBorders(size_t row) {
    double *rw = curr + row * width;

    MPI_Sendrecv(rw + (width - 2), 1, MPI_DOUBLE, rightN, TAG_LEFT_TO_RIGHT_ROW + (int)row,
                 rw, 1, MPI_DOUBLE, leftN, TAG_LEFT_TO_RIGHT_ROW + (int)row,
                 comm, MPI_STATUS_IGNORE);
}

void Field::sendFirstPass(size_t row) {
    if (rightN != NOBODY) {
        debug("s first pass");
        double bff[] = { bF[width - 2], cF[width - 2], fF[width - 2] };
        MPI_Send(bff, 3, MPI_DOUBLE, rightN, TAG_FIRST_PASS + (int)row, comm);
        debug("s first pass");
    }
}

void Field::receiveFirstPass(size_t row) {
    if (leftN != NOBODY) {
        debug("r first pass");
        double bff[3];
        MPI_Recv(bff, 3, MPI_DOUBLE, leftN, TAG_FIRST_PASS + (int)row, comm, MPI_STATUS_IGNORE);
        bF[0] = bff[0];
        cF[0] = bff[1];
        fF[0] = bff[2];
        debug("r first pass");
    }
}

void Field::sendSecondPass(size_t row) {
    if (leftN != NOBODY) {
        debug("s second pass");
        double *y = curr + row * width;
        MPI_Send(y + 1, 1, MPI_DOUBLE, leftN, TAG_SECOND_PASS + (int)row, comm);
        debug("s second pass");
    }
}

void Field::receiveSecondPass(size_t row) {
    if (rightN != NOBODY) {
        debug("r second pass");
        double *y = curr + row * width;
        MPI_Recv(y + width - 1, 1, MPI_DOUBLE, rightN, TAG_SECOND_PASS + (int)row, comm, MPI_STATUS_IGNORE);
        debug("r second pass");
    }
}

void Field::reduceMaxDelta(double &maxDelta) {
    MPI_Allreduce(&maxDelta, &maxDelta, 1, MPI_DOUBLE, MPI_MAX, rowComm);
}

void Field::reduceViews() {
    for (size_t index = 0, len = ftr.ViewCount(); index < len; ++index) {
        double value = view(index);
        double result = 0;
        MPI_Reduce(&value, &result, 1, MPI_DOUBLE, MPI_MAX, MASTER, comm);
        views[index] = result;
    }
}

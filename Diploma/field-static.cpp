//
//  Copyright © 2016 Nikolay Volosatov. All rights reserved.
//

#include "field-static.h"
#include "algo.h"
#include "balancing.h"
#include <cmath>

#include <sys/types.h>
#include <unistd.h>

void FieldStatic::init() {
    Field::init();
}

FieldStatic::~FieldStatic() {
    delete[] calculatingRows;
    delete[] nextCalculatingRows;

    delete[] sendBuff;
    delete[] receiveBuff;

    delete[] weights;

    delete[] nowBuckets;
    delete[] nextBuckets;
    
    MPI_Comm_free(&firstPassComm);
    MPI_Comm_free(&secondPassComm);
    MPI_Comm_free(&calculatingRowsComm);
}

void FieldStatic::finalize() {
    MPI_Barrier(MPI_COMM_WORLD);
}

void FieldStatic::calculateNBS() {
    Field::calculateNBS();

    fullHeight = algo::ftr().X2SplitCount();

    nowBuckets = new size_t[numProcs];
    nextBuckets = new size_t[numProcs];
    for (size_t i = 0; i < numProcs - 1; ++i) {
        nowBuckets[i] = height;
    }
    nowBuckets[numProcs - 1] = fullHeight - height * (numProcs - 1);

    shouldBalanceNext = false;
    shouldSendWeights = false;
    balancingCounter = (int)algo::ftr().TransposeBalanceIterationsInterval();

    if (bottomN == NOBODY) {
        height = fullHeight - height * (numProcs - 1);
    }

    height += (topN != NOBODY ? 1 : 0) + (bottomN != NOBODY ? 1 : 0);

    calculatingRows = new bool[width];
    nextCalculatingRows = new bool[width];

    sendBuff = new double[width * SEND_PACK_SIZE];
    receiveBuff = new double[width * SEND_PACK_SIZE];

    lastIterationsCount = lastWaitingCount = 0;

    MPI_Comm_dup(comm, &firstPassComm);
    MPI_Comm_dup(comm, &secondPassComm);
    MPI_Comm_dup(comm, &calculatingRowsComm);

    bundleSizeLimit = std::max(ceil((double)width / numProcs / 2), 15.0);
    printf("I'm %d(%d)\twith w:%zu\th:%zu\tbs:%zu.\tTop:%d\tbottom:%d\n",
           myId, ::getpid(), width, height, bundleSizeLimit, topN, bottomN);

    weights = new double[fullHeight];
    memset(weights, 0, fullHeight * sizeof(double));
}

#pragma mark - Logic

void FieldStatic::transpose(double **arr) {
    for (size_t index = 0, len = width * height; index < len; ++index) {
        size_t newIndex = (index % width) * height + index / width;
        buff[newIndex] = (*arr)[index];
    }
    std::swap(*arr, buff);
}

void FieldStatic::transpose() {
    transpose(transposed ? &curr : &prev);

    std::swap(hX, hY);
    std::swap(mySX, mySY);
    std::swap(topN, leftN);
    std::swap(bottomN, rightN);
    std::swap(width, height);
    transposed = transposed == false;
}

void FieldStatic::resetCalculatingRows() {
    for (size_t index = 0; index < height; ++index) {
        calculatingRows[index] = nextCalculatingRows[index] = true;
    }
}

size_t FieldStatic::firstPasses(size_t fromRow, bool first, bool async) {
    if (async && checkIncomingFirstPass(fromRow) == false) {
        return 0;
    }

    bool shouldProcess = recieveFirstPass(fromRow, first); // (crf + b + c + f) x [calculatingRows]

    if (shouldProcess == false) {
        // Sorry for that "workaround"
        return height * 2;
    }

    size_t row = fromRow;
    for (size_t bundleSize = 0; row < height && bundleSize < bundleSizeLimit; ++row) {
        if (calculatingRows[row] == false) {
            continue;
        }

        fillFactors(row, first);
        firstPass(row);
        ++bundleSize;
    }
    sendFirstPass(fromRow); // (crf + b + c + f) x [calculatingRows]

    return row;
}

size_t FieldStatic::secondPasses(size_t fromRow, bool first, bool async) {
    if (checkIncomingSecondPass(fromRow) == false) {
        if (async) {
            return 0;
        }
        else if (myCoord == 0) {
            ++lastWaitingCount;
        }
    }

    if (myCoord == 0) {
        ++lastIterationsCount;
    }

    recieveSecondPass(fromRow); // (nextCalculatingRows + y) x [prevCalculatingRows]

    size_t row = fromRow;
    for (size_t bundleSize = 0; row < height && bundleSize < bundleSizeLimit; ++row) {
        if (calculatingRows[row] == false) {
            nextCalculatingRows[row] = false;
            continue;
        }

        double delta = secondPass(row, first);

        nextCalculatingRows[row] = (rightN == NOBODY ? false : nextCalculatingRows[row]) || delta > epsilon;

        ++bundleSize;
    }
    sendSecondPass(fromRow); // (nextCalculatingRows + y) x [prevCalculatingRows]

    return row;
}

void FieldStatic::balanceBundleSize() {
    //printf("%zu\t%zu\n", lastWaitingCount, lastIterationsCount);
    if (lastWaitingCount > lastIterationsCount * algo::ftr().BalanceFactor()) {
        bundleSizeLimit = std::max(bundleSizeLimit - 1, algo::ftr().MinimumBundle());
    }
    else {
        bundleSizeLimit = std::min(bundleSizeLimit + 1, height / 2);
    }
    lastIterationsCount = 0;
    lastWaitingCount = 0;
}

size_t FieldStatic::solveRows() {
    size_t maxIterationsCount = 0;

    if (transposed) {
        --balancingCounter;
        if (balancingCounter == 0) {
            shouldSendWeights = true;
            balancingCounter = (int)algo::ftr().TransposeBalanceIterationsInterval();
        }

        resetCalculatingRows();
        bool first = true;

        if (myCoord == 0) {
            balanceBundleSize();
        }

        bool solving = true;
        while (solving) {
            size_t fromFirstPassRow = 0;
            size_t fromSecondPassRow = 0;

            while (fromSecondPassRow < height) {
                size_t nextSecondPassRow = 0;
                if (fromSecondPassRow < height && fromFirstPassRow > fromSecondPassRow) {
                    nextSecondPassRow = secondPasses(fromSecondPassRow, first, true);
                }

                size_t nextFirstPassRow = 0;
                if (fromFirstPassRow < height) {
                    nextFirstPassRow = firstPasses(fromFirstPassRow, first, true);
                }

                if (nextFirstPassRow == 0 && nextSecondPassRow == 0) {
                    if (fromFirstPassRow < height) {
                        nextFirstPassRow = firstPasses(fromFirstPassRow, first, false);
                    }
                    else {
                        nextSecondPassRow = secondPasses(fromSecondPassRow, first, false);
                    }
                }

                if (nextFirstPassRow > 0) {
                    if (nextFirstPassRow == height * 2) {
                        solving = false;
                        break;
                    }
                    fromFirstPassRow = nextFirstPassRow;
                }
                if (nextSecondPassRow > 0) {
                    fromSecondPassRow = nextSecondPassRow;
                }
            }

            if (solving == false) {
                break;
            }

            first = false;
            std::swap(nextCalculatingRows, calculatingRows);
            ++maxIterationsCount;
            if (maxIterationsCount >= MAX_ITTERATIONS_COUNT) {
                break;
            }

            if (leftN == NOBODY) {
                solving = false;
                for (size_t row = 0; row < height; ++row) {
                    if (calculatingRows[row]) {
                        solving = true;
                        break;
                    }
                }

                if (solving == false) {
                    sendDoneAsFirstPass();
                }
            }
        }
    }
    else {
        if (shouldBalanceNext) {
            balance();
            shouldBalanceNext = false;
        }

        size_t firstRealRow = (topN == NOBODY ? 0 : 1);
        size_t lastRealRow = height - firstRealRow - (bottomN ? 0 : 1);
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
            maxIterationsCount = std::max(maxIterationsCount, iterationsCount);

            if (firstRealRow <= row && row <= lastRealRow) {
                weights[mySY + row - firstRealRow] =
                        weights[mySY + row - firstRealRow] * algo::ftr().TransposeBalanceFactor()
                        + iterationsCount * (1.0 - algo::ftr().TransposeBalanceTimeFactor())
                        + (picosecFromStart() - startTime) * 1e-12 * algo::ftr().TransposeBalanceTimeFactor();
            }
        }
    }
    
    return maxIterationsCount;
}

#pragma mark - MPI

void FieldStatic::sendFirstPass(size_t fromRow) {
    // (crf + b + c + f) x [calculatingRows]
    if (rightN != NOBODY) {
        double *sBuff = sendBuff + fromRow * SEND_PACK_SIZE;

        int sSize = 0;
        sBuff[sSize++] = (shouldSendWeights ? -1 : 1) * (int)bundleSizeLimit;

        if (shouldSendWeights) {
            for (size_t i = 0; i < mySX + width - (rightN == NOBODY ? 0 : 1) - (leftN == NOBODY ? 0 : 1); ++i) {
                sBuff[sSize++] = weights[i];
            }
            shouldSendWeights = false;
        }

        for (size_t row = fromRow, bundleSize = 0; row < height && bundleSize < bundleSizeLimit; ++row) {
            if (calculatingRows[row] == false) {
                continue;
            }

            size_t index = row * width + width - 2;
            sBuff[sSize++] = row;
            sBuff[sSize++] = mbF[index];
            sBuff[sSize++] = mcF[index];
            sBuff[sSize++] = mfF[index];
            ++bundleSize;
        }

        MPI_Request request;
        MPI_Isend(sBuff, sSize, MPI_DOUBLE, rightN, (int)fromRow, firstPassComm, &request);
    }
}

void FieldStatic::sendDoneAsFirstPass() {
    if (rightN != NOBODY) {
        MPI_Request request;
        MPI_Isend(NULL, 0, MPI_DOUBLE, rightN, 0, firstPassComm, &request);
    }
}

bool FieldStatic::checkIncomingFirstPass(size_t fromRow) {
    if (leftN == NOBODY) {
        return true;
    }

    int flag;
    MPI_Iprobe(leftN, (int)fromRow, firstPassComm, &flag, MPI_STATUS_IGNORE);
    return flag;
}

bool FieldStatic::recieveFirstPass(size_t fromRow, bool first) {
    // (crf + b + c + f) x [calculatingRows]
    if (leftN != NOBODY) {
        MPI_Status status;
        MPI_Probe(leftN, (int)fromRow, firstPassComm, &status);
        int sSize;
        MPI_Get_count(&status, MPI_DOUBLE, &sSize);

        MPI_Request request;
        MPI_Irecv(receiveBuff, sSize, MPI_DOUBLE, leftN, (int)fromRow, firstPassComm, &request);
        if (sSize == 0) {
            MPI_Cancel(&request);
            sendDoneAsFirstPass();
            return false;
        } else {
            MPI_Wait(&request, MPI_STATUS_IGNORE);
        }

        size_t idxBuffer = 0;
        double bsFlag = receiveBuff[idxBuffer++];
        bool receiveWeights = bsFlag < 0;
        bundleSizeLimit = (receiveWeights ? -1 : 1) * bsFlag;

        if (receiveWeights) {
            for (size_t i = 0; i < mySX; ++i) {
                weights[i] = receiveBuff[idxBuffer++];
            }
            shouldSendWeights = true;

            if (rightN == NOBODY) {
                partitionAndCheck();
            }
        }

        size_t lastRow = fromRow;
        for (; idxBuffer < sSize;) {
            size_t row = receiveBuff[idxBuffer++];

            while (lastRow < row) {
                calculatingRows[lastRow++] = false;
            }
            calculatingRows[lastRow++] = true;

            size_t index = row * width;
            mbF[index] = receiveBuff[idxBuffer++];
            mcF[index] = receiveBuff[idxBuffer++];
            mfF[index] = receiveBuff[idxBuffer++];
        }

        if ((sSize - 1) / 4 < bundleSizeLimit) {
            while (lastRow < height) {
                calculatingRows[lastRow++] = false;
            }
        }
    }

    return true;
}

void FieldStatic::partitionAndCheck() {
    /*debug() << "NB ";
    for (size_t i = 0; i < numProcs; ++i) {
        debug(0) << nowBuckets[i] << " ";
    }
    (debug(0) << "\n").flush();*/

    auto buckets = balancing::fastPartition(weights, fullHeight, nowBuckets, numProcs);
    size_t deltaSum = 0;
    for (size_t i = 0; i < numProcs; ++i) {
        nextBuckets[i] = (size_t)buckets[i];
        deltaSum += std::abs((long)nextBuckets[i] - (long)nowBuckets);
    }
    // TODO: Add config parameter
    shouldBalanceNext = deltaSum > 3;
}

void FieldStatic::sendSecondPass(size_t fromRow) {
    // (nextCalculatingRows + y) x [prevCalculatingRows]
    if (leftN != NOBODY) {
        double *sBuff = sendBuff + fromRow * SEND_PACK_SIZE;
        int sSize = 0;

        sBuff[sSize++] = shouldBalanceNext ? 1 : 0;
        if (shouldBalanceNext) {
            for (size_t i = 0; i < numProcs; ++i) {
                sBuff[sSize++] = nextBuckets[i];
            }
        }

        for (size_t row = fromRow, bundleSize = 0; row < height && bundleSize < bundleSizeLimit; ++row) {
            if (calculatingRows[row] == false) {
                continue;
            }
            sBuff[sSize++] = (nextCalculatingRows[row] ? 1 : -1) * curr[row * width + 1];

            ++bundleSize;
        }

        MPI_Request request;
        MPI_Isend(sBuff, sSize, MPI_DOUBLE, leftN, (int)fromRow, secondPassComm, &request);
    }
}

bool FieldStatic::checkIncomingSecondPass(size_t fromRow) {
    if (rightN == NOBODY) {
        return true;
    }

    int flag;
    MPI_Iprobe(rightN, (int)fromRow, secondPassComm, &flag, MPI_STATUS_IGNORE);
    return flag;
}

void FieldStatic::recieveSecondPass(size_t fromRow) {
    // (nextCalculatingRows + y) x [prevCalculatingRows]
    if (rightN != NOBODY) {
        MPI_Status status;
        MPI_Probe(rightN, (int)fromRow, secondPassComm, &status);
        int sSize;
        MPI_Get_count(&status, MPI_DOUBLE, &sSize);

        MPI_Recv(receiveBuff, sSize, MPI_DOUBLE, rightN, (int)fromRow, secondPassComm, MPI_STATUS_IGNORE);

        sSize = 0;
        bool shouldReceiveBuckets = receiveBuff[sSize++] > 0;
        if (shouldReceiveBuckets) {
            for (size_t i = 0; i < numProcs; ++i) {
                nextBuckets[i] = receiveBuff[sSize++];
            }
            shouldBalanceNext = true;
        }

        for (size_t row = fromRow, bundleSize = 0; row < height && bundleSize < bundleSizeLimit; ++row) {
            if (calculatingRows[row] == false) {
                continue;
            }

            double value = receiveBuff[sSize++];
            nextCalculatingRows[row] = value > 0;
            curr[(row + 1) * width - 1] = nextCalculatingRows[row] ? value : -value;

            ++bundleSize;
        }
    }
}

#pragma mark - Balancing

bool FieldStatic::balanceNeeded() {
    // Balance in solving
    return false;
}

void FieldStatic::balance() {
    debug() << "Balancing " << t << "\n";
    debug(0).flush();
}

#pragma mark - Print

void FieldStatic::printConsole() {
    if (algo::ftr().EnableConsole()) {
        double viewValue = view(algo::ftr().DebugView());

        if (fabs(viewValue - NOTHING) > __DBL_EPSILON__) {
            printf("Field[%d] (itrs: %zu, bsL %zu, time: %.5f)\tview: %.7f\n",
                   myId, lastIterrationsCount, bundleSizeLimit, t, viewValue);
        }
    }
}

//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#include "field.h"
#include "algo.h"
#include <cmath>
#include <sys/types.h>
#include <unistd.h>
#include <fstream>
#include <limits>
#include <ios>
#include <iostream>

#ifdef DEBUG
#define DEBUG_PRINT
//#define DEBUG_WAIT
#endif

template<typename Ch, typename Traits = std::char_traits<Ch> >
struct basic_nullbuf : std::basic_streambuf<Ch, Traits> {
    typedef std::basic_streambuf<Ch, Traits> base_type;
    typedef typename base_type::int_type int_type;
    typedef typename base_type::traits_type traits_type;

    virtual int_type overflow(int_type c) {
        return traits_type::not_eof(c);
    }
};

std::ostream &debugStream(int myId) {
#ifdef DEBUG_PRINT
    char buf[255] = {0};
    sprintf(buf, "log.%d.txt", myId);
    static std::ofstream f(buf);
    return f;
#else
    static basic_nullbuf<char> nullbuf;
    static std::ostream cnull(&nullbuf);
    return cnull;
#endif
}

unsigned long long Field::picosecFromStart() {
    return std::chrono::duration<unsigned long long, std::pico>(bx_clock_t::now() - startSyncTime).count();
}

std::ostream &Field::debug(bool info) {
    std::ostream &out = debugStream(myId);
    if (info == false) {
        return out;
    }

    out.precision(15);
    out << "N[" << myId << " " << std::fixed << picosecFromStart() * 1e-12 << "] ";
    out.unsetf(std::ios_base::floatfield);
    out.precision(7);
    return out;
}

void Field::calculateNBS() {
    MPI_Barrier(MPI_COMM_WORLD);
    startSyncTime = bx_clock_t::now();

    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    int dims[] = { numProcs };
    int wrap[] = { 0 };
    MPI_Cart_create(MPI_COMM_WORLD, 1, dims, wrap, 1, &comm);

    MPI_Comm_rank(comm, &myId);
    MPI_Cart_coords(comm, myId, 1, &myCoord);
    MPI_Cart_shift(comm, 0, 1, &topN, &bottomN);
    rightN = leftN = NOBODY;

    height /= numProcs;

    mySY = height * myCoord;
    mySX = 0;

#ifdef DEBUG_WAIT
    printf("I'm %d(%d)\n", myId, ::getpid());
    int waiter = myId;
    while (waiter == WAITER) sleep(5);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void Field::finalize() {
}

void Field::reduceViews() {
    for (size_t index = 0, len = algo::ftr().ViewCount(); index < len; ++index) {
        double value = view(index);
        double result = 0;
        MPI_Reduce(&value, &result, 1, MPI_DOUBLE, MPI_MAX, MASTER, comm);
        views[index] = result;
    }
}

void Field::printMatrix() {
    if (algo::ftr().EnableMatrix()) {
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
                delete mfout;
                mfout = NULL;
            }
            MPI_Barrier(comm);
        }
    }
}

#pragma mark - Balancing

void Field::syncWeights() {
    
}

bool Field::balanceNeeded() {
    return false;
}

void Field::balance() {

}

bool Field::isBucketsMaster() {
    return false;
}

size_t Field::weightsSize() {
    return false;
}

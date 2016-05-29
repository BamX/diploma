//
//  Copyright (c) 2016 Nikolay Volosatov. All rights reserved.
//

#ifndef field_transpose_h
#define field_transpose_h

#include "field.h"

class FieldTranspose : public Field {

    MPI_Comm balanceComm;

    int balancingCounter;

    size_t mySYT;
    size_t *hBuckets, *vBuckets;
    std::vector<int> nextBuckets, nextBucketsT;

    MPI_Datatype *sendtypes, *recvtypes;
    int *sendcounts, *senddispls, *recvcounts, *recvdispls, *gathercounts, *gatherdispls;

    void transpose(double *arr);
    void transpose() override;

    void calculateNBS() override;

    size_t solveRows() override;

    void printConsole() override;
    void printMatrix() override;

#pragma mark - Balancing MPI

    double *weights, *weightsT;

    size_t heightCapacity;
    void resize(size_t newHeight);

    void createVType(size_t width, size_t height, size_t bWidth, MPI_Datatype *type);
    void createHType(size_t width, size_t height, size_t bWidth, MPI_Datatype *type);

    void syncWeights() override;
    bool balanceNeeded() override;
    void balance() override;

public:
    FieldTranspose(const char *filename);
    ~FieldTranspose();

    void init() override;
    double view(double x1, double x2) override;
};

#endif /* field_transpose_h */

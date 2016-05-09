//
//  Copyright (c) 2016 Nikolay Volosatov. All rights reserved.
//

#ifndef field_transpose_h
#define field_transpose_h

#include "field.h"

class FieldTranspose : public Field {

    MPI_Datatype mpiAllType;
    MPI_Comm balanceComm;

    void transpose(double *arr);
    void transpose() override;

    void calculateNBS() override;

    size_t solveRows() override;

    void printConsole() override;
    void printMatrix() override;

#pragma mark - Balancing MPI

    void syncWeights() override;
    bool balanceNeeded() override;
    void balance() override;

public:
    ~FieldTranspose();

    void init() override;
    double view(double x1, double x2) override;
};

#endif /* field_transpose_h */

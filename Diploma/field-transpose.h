//
//  Copyright (c) 2016 Nikolay Volosatov. All rights reserved.
//

#ifndef field_transpose_h
#define field_transpose_h

#include "field.h"

class FieldTranspose : public Field {

    MPI_Datatype mpiAllType;

    void transpose(double *arr);
    void transpose() override;

    size_t solveRows() override;

    void printConsole() override;

public:
    ~FieldTranspose();

    void init() override;
};

#endif /* field_transpose_h */

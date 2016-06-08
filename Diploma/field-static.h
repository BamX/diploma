//
//  Copyright © 2016 Nikolay Volosatov. All rights reserved.
//

#ifndef field_static_h
#define field_static_h

#include "field.h"

class FieldStatic : public Field {
    size_t bundleSizeLimit;
    size_t lastWaitingCount, lastIterationsCount;

    bool *calculatingRows, *nextCalculatingRows;
    double *sendBuff, *receiveBuff;

    void calculateNBS() override;
    void resetCalculatingRows();

    size_t firstPasses(size_t fromRow, bool first, bool async);
    size_t secondPasses(size_t fromRow, bool first, bool async);

    void transpose(double **arr);
    void transpose() override;

    size_t solveRows() override;

    void sendRecieveCalculatingRows();
    void balanceBundleSize();

    MPI_Comm firstPassComm, secondPassComm, calculatingRowsComm;

    void sendFirstPass(size_t fromRow);
    void sendDoneAsFirstPass();
    bool checkIncomingFirstPass(size_t fromRow);
    bool recieveFirstPass(size_t fromRow, bool first);
    
    void sendSecondPass(size_t fromRow);
    bool checkIncomingSecondPass(size_t fromRow);
    void recieveSecondPass(size_t fromRow);

    void printConsole() override;
    
public:
    ~FieldStatic();

    void init() override;
    void finalize() override;
};

#endif /* field_static_h */

//
//  Copyright © 2016 Nikolay Volosatov. All rights reserved.
//

#ifndef field_static_h
#define field_static_h

#include "field.h"

class FieldStatic : public Field {
    size_t bundleSizeLimit;
    size_t lastWaitingCount, lastIterationsCount;

    void resetCalculatingRows();

    size_t firstPasses(size_t fromRow, bool first, bool async);
    size_t secondPasses(size_t fromRow, bool first, bool async);

    size_t solveRows() override;

    void sendRecieveCalculatingRows();
    void balanceBundleSize();

    MPI_Comm firstPassComm, secondPassComm, calculatingRowsComm;

    void sendFistPass(size_t fromRow);
    bool checkIncomingFirstPass(size_t fromRow);
    void recieveFirstPass(size_t fromRow, bool first);
    void sendSecondPass(size_t fromRow);
    bool checkIncomingSecondPass(size_t fromRow);
    void recieveSecondPass(size_t fromRow);

    void printConsole() override;
    
public:
    FieldStatic();
    ~FieldStatic();
};

#endif /* field_static_h */

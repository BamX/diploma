//
//  Copyright (c) 2016 Nikolay Volosatov. All rights reserved.
//

#ifndef algo_h
#define algo_h

#include "factors.h"

namespace algo {

    Factors &ftr();

    void fillFactors(double *rw, double *brw, size_t size,
                     double *aF, double *bF, double *cF, double *fF,
                     double t, double hX, double dT,
                     bool leftBorder, bool rightBorder);

    void firstPass(size_t size, double *aF, double *bF, double *cF, double *fF);
    
    void secondPass(double *rw, double *brw, size_t size,
                    double *bF, double *cF, double *fF,
                    bool rightBorder, double *maxDelta);

}

#endif /* algo_h */

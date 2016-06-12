#include "algo.h"
#include "factors.h"
#include <cmath>
#include <algorithm>

namespace algo {

    static Factors factors;

    Factors &ftr() {
        return factors;
    }

    void fillFactors(double *rw, double *brw, size_t size,
                     double *aF, double *bF, double *cF, double *fF,
                     double t, double hX, double dT,
                     bool leftBorder, bool rightBorder)
    {
        double lm0 = factors.lambda(brw[0]), lmh = factors.lambda(brw[1]);
        if (leftBorder) {
            aF[0] = 0;
            cF[0] = dT * (lm0 + lmh) + hX * hX * factors.ro(brw[0]) * factors.cEf(brw[0]);
            bF[0] = -dT * (lm0 + lmh);
            fF[0] = hX * hX * factors.ro(brw[0]) * factors.cEf(brw[0]) * rw[0];
        }

        if (rightBorder) {
            double TPrev = brw[size - 1];
            double TPrev4 = TPrev * TPrev * TPrev * TPrev;
            double lmXX = factors.lambda(brw[size - 1]), lmXXm1 = factors.lambda(brw[size - 2]);

            aF[size - 1] = -dT * (lmXXm1 + lmXX);
            cF[size - 1] = dT * (lmXXm1 + lmXX)
                    + hX * hX * factors.ro(brw[size - 1]) * factors.cEf(brw[size - 1])
                    + 2 * hX * dT * factors.alpha(t);
            bF[size - 1] = 0;
            fF[size - 1] = hX * hX * factors.ro(brw[size - 1]) * factors.cEf(brw[size - 1]) * rw[size - 1]
                    - 2 * hX * dT * factors.sigma(t) * (TPrev4 - factors.TEnv4())
                    + 2 * hX * dT * factors.alpha(t) * factors.TEnv();
        }

        double lmXm1 = lm0, lmX = lmh, lmXp1;
        double mhh2rocdT;
        for (size_t index = 1; index < size - 1; ++index) {
            lmXp1 = factors.lambda(brw[index + 1]);
            mhh2rocdT = - 2 * hX * hX * factors.ro(brw[index]) * factors.cEf(brw[index]) / dT;

            aF[index] = lmX + lmXm1;
            bF[index] = lmXp1 + lmX;
            cF[index] = -(lmXp1 + 2 * lmX + lmXm1) + mhh2rocdT;
            fF[index] = mhh2rocdT * rw[index];
            
            lmXm1 = lmX;
            lmX = lmXp1;
        }
    }

    void firstPass(size_t size, double *aF, double *bF, double *cF, double *fF, bool rightBorder) {
        double m;
        for (size_t i = 1, len = size - (rightBorder ? 0 : 1); i < len; ++i) {
            m = aF[i] / cF[i - 1];
            cF[i] -= m * bF[i - 1];
            fF[i] -= m * fF[i - 1];
        }
    }

    void secondPass(double *rw, double *brw, size_t size,
                    double *bF, double *cF, double *fF,
                    bool rightBorder, double *maxDelta) {
        *maxDelta = 0;

        double newValue = 0;
        if (rightBorder) {
            newValue = fF[size - 1] / cF[size - 1];
            *maxDelta = fabs(newValue - rw[size - 1]);
            brw[size - 1] = newValue;
        }

        for (long i = size - 2; i >= 0; --i) {
            newValue = (fF[i] - bF[i] * brw[i + 1]) / cF[i];

            double newDelta = fabs(newValue - rw[i]);
            *maxDelta = std::max(*maxDelta, newDelta);
            brw[i] = newValue;
        }
    }

}
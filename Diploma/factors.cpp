//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#include "factors.h"
#include <cmath>

namespace ftr {
    static double const moveVelocity = 0.75 / 60; // м/с

    static double const TLik = 1738; // К
    static double const TSol = 1679; // К

    static double const L = 272; // кДж/кг
    static double const cLik = 710; // Дж/(кг * К)

    double const X1 = 0.15;  // м
    double const X2 = 0.125; // м
    double const totalTime = (0.4 + 0.4 + 0.47 + 0.95 + 1.51 + 18.97) / moveVelocity; // с

    double const TStart = 1768; // К
    double const TEnv = 303; // К
    double const TEnv4 = TEnv * TEnv * TEnv * TEnv; // К

    inline double Li(unsigned short i, double x) {
        switch (i) {
            case 1:
                return 44076 - 85622 * x * x + 50357 * x;
            case 2:
                return 5163.2 - 74009 * x * x + 70232 * x;
            case 3:
                return 2622.3 - 92590 * x * x + 80523 * x;
            case 4:
                return 14775 - 154544 * x * x + 142489 * x;
        }
        return 0.0;
    }

    inline double Ti(unsigned short i) {
        switch (i) {
            case 1:
                return 1000;
            case 2:
                return 1033;
            case 3:
                return 923;
            case 4:
                return 1033;
        }
        return 0.0;
    }

    inline double dTi(unsigned short i) {
        switch (i) {
            case 1:
                return 70;
            case 2:
                return 350;
            case 3:
                return 1100;
            case 4:
                return 170;
        }
        return 0.0;
    }

    inline double cSol(double T) {
        double result = 469 + 0.16 * (T - 323);
        for (unsigned short i = 1; i <= 4; ++i) {
            double tC = (Ti(i) - T) / dTi(i);
            result += 4.5141 * Li(i, T) / dTi(i) * exp(-16 * tC * tC);
        }
        return result;
    }

    inline double sigm(double T) {
        double cLikS = 15.463359 - 0.124528e-1 * T + 0.216279e-5 * T * T;
        double cSolS = -11.0388 + 0.278656e-1 * T - 0.120163e-4 * T * T;
        return (cLikS - 0.7) / (cLikS - cSolS);
    }

    double cEf(double T, double dT) {
        if (T >= ftr::TLik) {
            return ftr::cLik;
        }
        else if (T > ftr::TSol) {
            return ftr::cSol(T) - ftr::L * (ftr::sigm(T + dT) - ftr::sigm(T)) / dT;
        }
        else {
            return ftr::cSol(T);
        }
    }

    double alpha(double t) {
        double x = t * moveVelocity;

        if (x <= 0.4) {
            return 2100;
        }
        else if (x <= 0.4 + 0.4) {
            return 60;
        }
        else if (x <= 0.4 + 0.4 + 0.47) {
            return 850;
        }
        else if (x <= 0.4 + 0.4 + 0.47 + 0.95) {
            return 120;
        }
        else if (x <= 0.4 + 0.4 + 0.47 + 0.95 + 1.51) {
            return 40;
        }
        else {
            return 25;
        }
    }

    double sigma(double t) {
        double x = t * moveVelocity;

        if (x <= 0.4 + 0.4 + 0.47 + 0.95 + 1.51) {
            return 0;
        }
        else {
            return 3.2e-8;
        }
    }
}

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
    static double const x = 0.7;

    static double const totalLength = (0.4 + 0.4 + 0.47 + 0.95 + 1.51 + 18.97); // с

    static double Temps[] = { 273, 373, 473, 573, 673, 773, 873,
        973, 1073, 1173, 1273, 1373, 1473, 1679, 1682, 1800
    };
    static double Lambds[] = {
        52.56057, 51.35258, 49.16971, 46.22939, 42.74907, 38.94618, 35.03819,
        31.24254, 24.06517, 25.37233, 26.95363, 28.32515, 29.40302, 31.62343,
        28.0, 28.0
    };
    static double Ros[] = {
        7885.884, 7845.138, 7804.392, 7763.646, 7722.9, 7682.901, 7647.512, 7621.141,
        7631.934, 7572.359, 7512.151, 7453.459, 7398.322, 7190.562, 7000.0, 7000.0,
    };

    static double Ti[] = { 1000, 1033, 923, 1033 };
    static double dTi[] = { 70, 350, 1100, 170 };

    inline double alphaForT(double T, size_t &index) {
        index = 0;
        while (Temps[index] < T) ++index;
        return (T - Temps[index - 1]) / (Temps[index] - Temps[index - 1]);
    }

    inline double Lambda(double T) {
        size_t index = 0;
        double alpha = alphaForT(T, index);
        return Lambds[index - 1] + alpha * (Lambds[index] - Lambds[index - 1]);
    }

    inline double Ro(double T) {
        size_t index = 0;
        double alpha = alphaForT(T, index);
        return Ros[index - 1] + alpha * (Ros[index] - Ros[index - 1]);
    }

    inline double Li(unsigned short i, double x) {
        switch (i) {
            case 0:
                return 44076 - 85622 * x * x + 50357 * x;
            case 1:
                return 5163.2 - 74009 * x * x + 70232 * x;
            case 2:
                return 2622.3 - 92590 * x * x + 80523 * x;
            case 3:
                return 14775 - 154544 * x * x + 142489 * x;
        }
        return 0.0;
    }

    inline double cSol(double T) {
        double result = 469 + 0.16 * (T - 323);
        for (unsigned short i = 0; i < 4; ++i) {
            double tC = (Ti[i] - T) / dTi[i];
            result += 4.5141 * Li(i, x) / dTi[i] * exp(-16 * tC * tC);
        }
        return result;
    }

    inline double sigm(double T) {
        double cLikS = 15.463359 - 0.124528e-1 * T + 0.216279e-5 * T * T;
        double cSolS = -11.0388 + 0.278656e-1 * T - 0.120163e-4 * T * T;
        return (cLikS - 0.7) / (cLikS - cSolS);
    }

    inline double ksiFunc(double T) { // k = 0.7
        double x2 = x * x;
        double z = (1.8691e6 - 2843.51*x + x*x);
        return (1.31914e9 - 1.51221e6*x + 444.52*x2)/(z * z);
    }
}

Factors::Factors() {
    _x1 = config.value("X1");
    _x2 = config.value("X2");
    _totalTime = ftr::totalLength / config.value("Speed");
    _x1SplitCount = config.value("X1SplitCount");
    _x2SplitCount = config.value("X2SplitCount");
    _timeSplitCount = config.value("TimeSplitCount");
    _epsilon = config.value("Epsilon");

    _TStart = config.value("InitT");
    _TEnv = config.value("EnvT");
    _TEnv4 = _TEnv * _TEnv * _TEnv * _TEnv;

    _enableConsole = config.value("EnableConsole") > 0;
    _enablePlot = config.value("EnablePlot") > 0;
    _enableMatrix = config.value("EnableMatrix") > 0;

    _viewCount = config.value("ViewCount");
    _debugView = config.value("DebugView");
    _matrixFramesCount = config.value("MatrixFramesCount");

    for (size_t index = 0; index < _viewCount; ++index) {
        char buff[10];
        snprintf(buff, 10, "View%zuX1", index);
        _x1View.push_back(config.value(buff));

        snprintf(buff, 10, "View%zuX2", index);
        _x2View.push_back(config.value(buff));
    }
}

double Factors::cEf(double T) const {
    if (T >= ftr::TLik) {
        return ftr::cLik;
    }
    else if (T > ftr::TSol) {
        return ftr::cSol(T) - ftr::L * ftr::ksiFunc(T);
    }
    else {
        return ftr::cSol(T);
    }
}

double Factors::alpha(double t) const {
    double x = t * ftr::moveVelocity;

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

double Factors::sigma(double t) const {
    double x = t * ftr::moveVelocity;

    if (x <= 0.4 + 0.4 + 0.47 + 0.95 + 1.51) {
        return 0;
    }
    else {
        return 3.2e-8;
    }
}

double Factors::lambda(double T) const {
    return ftr::Lambda(T);
}

double Factors::ro(double T) const {
    return ftr::Ro(T);
}

double Factors::X1() const {
    return _x1;
}

double Factors::X2() const {
    return _x2;
}

double Factors::totalTime() const {
    return _totalTime;
}

double Factors::X1SplitCount() const {
    return _x1SplitCount;
}

double Factors::X2SplitCount() const {
    return _x2SplitCount;
}

double Factors::TimeSplitCount() const {
    return _timeSplitCount;
}

double Factors::Epsilon() const {
    return _epsilon;
}

double Factors::TStart() const {
    return _TStart;
}

double Factors::TEnv() const {
    return _TEnv;
}

double Factors::TEnv4() const {
    return _TEnv4;
}

bool Factors::EnableConsole() const {
    return _enableConsole;
}

bool Factors::EnablePlot() const {
    return _enablePlot;
}

bool Factors::EnableMatrix() const {
    return _enableMatrix;
}

size_t Factors::ViewCount() const {
    return _viewCount;
}

size_t Factors::DebugView() const {
    return _debugView;
}

size_t Factors::MatrixFramesCount() const {
    return _matrixFramesCount;
}

double Factors::X1View(size_t index) const {
    return _x1View[index];
}

double Factors::X2View(size_t index) const {
    return _x2View[index];
}

//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#include "factors.h"
#include <cmath>

namespace ftr {
    static float const moveVelocity = 0.75 / 60; // м/с

    static float const TLik = 1738; // К
    static float const TSol = 1679; // К
    static float const dT = 0.000001; // K ???

    static float const L = 272; // кДж/кг
    static float const cLik = 710; // Дж/(кг * К)
    static float const x = 0.7;

    static float const totalLength = (0.4 + 0.4 + 0.47 + 0.95 + 1.51 + 18.97); // с

    static float Temps[] = { 273, 373, 473, 573, 673, 773, 873,
        973, 1073, 1173, 1273, 1373, 1473, 1679, 1682, 1800
    };
    static float Lambds[] = {
        52.56057, 51.35258, 49.16971, 46.22939, 42.74907, 38.94618, 35.03819,
        31.24254, 24.06517, 25.37233, 26.95363, 28.32515, 29.40302, 31.62343,
        28.0, 28.0
    };
    static float Ros[] = {
        7885.884, 7845.138, 7804.392, 7763.646, 7722.9, 7682.901, 7647.512, 7621.141,
        7631.934, 7572.359, 7512.151, 7453.459, 7398.322, 7190.562, 7000.0, 7000.0,
    };

    static float Ti[] = { 1000, 1033, 923, 1033 };
    static float dTi[] = { 70, 350, 1100, 170 };

    inline float alphaForT(float T, size_t &index) {
        index = 0;
        while (Temps[index] < T) ++index;
        return (T - Temps[index - 1]) / (Temps[index] - Temps[index - 1]);
    }

    inline float Lambda(float T) {
        size_t index = 0;
        float alpha = alphaForT(T, index);
        return Lambds[index - 1] + alpha * (Lambds[index] - Lambds[index - 1]);
    }

    inline float Ro(float T) {
        size_t index = 0;
        float alpha = alphaForT(T, index);
        return Ros[index - 1] + alpha * (Ros[index] - Ros[index - 1]);
    }

    inline float Li(unsigned short i, float x) {
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

    inline float cSol(float T) {
        float result = 469 + 0.16 * (T - 323);
        for (unsigned short i = 0; i < 4; ++i) {
            float tC = (Ti[i] - T) / dTi[i];
            result += 4.5141 * Li(i, x) / dTi[i] * exp(-16 * tC * tC);
        }
        return result;
    }

    inline float sigm(float T) {
        float cLikS = 15.463359 - 0.124528e-1 * T + 0.216279e-5 * T * T;
        float cSolS = -11.0388 + 0.278656e-1 * T - 0.120163e-4 * T * T;
        return (cLikS - 0.7) / (cLikS - cSolS);
    }
}

float Factors::cEf(float T) const {
    if (T >= ftr::TLik) {
        return ftr::cLik;
    }
    else if (T > ftr::TSol) {
        return ftr::cSol(T) - ftr::L * (ftr::sigm(T + ftr::dT) - ftr::sigm(T)) / ftr::dT;
    }
    else {
        return ftr::cSol(T);
    }
}

float Factors::alpha(float t) const {
    float x = t * ftr::moveVelocity;

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

float Factors::sigma(float t) const {
    float x = t * ftr::moveVelocity;

    if (x <= 0.4 + 0.4 + 0.47 + 0.95 + 1.51) {
        return 0;
    }
    else {
        return 3.2e-8;
    }
}

float Factors::lambda(float T) const {
    return ftr::Lambda(T);
}

float Factors::ro(float T) const {
    return ftr::Ro(T);
}

float Factors::X1() const {
    return config.value("X1");
}

float Factors::X2() const {
    return config.value("X2");
}

float Factors::totalTime() const {
    return ftr::totalLength / config.value("Speed");
}

float Factors::X1SplitCount() const {
    return config.value("X1SplitCount");
}

float Factors::X2SplitCount() const {
    return config.value("X2SplitCount");
}

float Factors::TimeSplitCount() const {
    return config.value("TimeSplitCount");
}

float Factors::Epsilon() const {
    return config.value("Epsilon");
}

float Factors::TStart() const {
    return config.value("InitT");
}

float Factors::TEnv() const {
    return config.value("EnvT");
}

float Factors::TEnv4() const {
    float value = TEnv();
    value *= value;
    value *= value;
    return value;
}

bool Factors::EnableConsole() const {
    return config.value("EnableConsole") > 0;
}

bool Factors::EnablePlot() const {
    return config.value("EnablePlot") > 0;
}

bool Factors::EnableMatrix() const {
    return config.value("EnableMatrix") > 0;
}

size_t Factors::ViewCount() const {
    return config.value("ViewCount");
}

size_t Factors::DebugView() const {
    return config.value("DebugView");
}

size_t Factors::MatrixFramesCount() const {
    return config.value("MatrixFramesCount");
}

float Factors::X1View(size_t index) const {
    char buff[10];
    snprintf(buff, 10, "View%luX1", index);
    std::string key = buff;

    return config.value(key);
}

float Factors::X2View(size_t index) const {
    char buff[10];
    snprintf(buff, 10, "View%luX2", index);
    std::string key = buff;

    return config.value(key);
}

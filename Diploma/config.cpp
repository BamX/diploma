//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#include "config.h"

#include <fstream>

Config::Config() {
    std::ifstream fin("config.ini");

    std::string name;
    double value;

    while (fin.good()) {
        fin >> name >> value;
        if (name.length() > 0) {
            values[name] = value;
        }
    }

    fin.close();
}

double Config::value(std::string name) const {
    return values.at(name);
}

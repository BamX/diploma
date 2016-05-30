//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#include "config.h"

#include <fstream>

Config::Config(const char *filename) {
    std::ifstream fin(filename ?: "config.ini");

    std::string name;
    std::string value;

    while (fin.good()) {
        fin >> name >> value;
        if (name.length() > 0 && name[0] != '#') {
            values[name] = value;
        }
    }

    fin.close();
}

double Config::value(std::string name) const {
    return atof(values.at(name).data());
}

std::string Config::str_value(std::string name) const {
    return values.at(name);
}

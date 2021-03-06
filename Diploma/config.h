//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#ifndef __Diploma__config__
#define __Diploma__config__

#include <map>
#include <string>

class Config {
    std::map<std::string, std::string> values;
public:

    Config(const char *filename);
    double value(std::string name) const;
    std::string str_value(std::string name) const;
};

#endif /* defined(__Diploma__config__) */

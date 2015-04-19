//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#include <stdio.h>

#include "field.h"

int main(int argc, const char * argv[]) {
    Field field(60, 60, 1.2);
    field.fillInitial();
    field.print();

    size_t limit = 50;
    while (field.done() == false && limit --> 0) {
        field.solve();
        field.print();
    }

    return 0;
}

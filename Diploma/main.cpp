//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#include <stdio.h>

#include "field.h"

int main(int argc, const char * argv[]) {
    Field field(100, 100, 200000);
    field.test();
    field.enableFileOutput();
    
    field.fillInitial();
    field.print();

    size_t limit = 20;//INT_MAX;
    while (field.done() == false && limit --> 0) {
        field.solve();
        field.print();
    }

    return 0;
}

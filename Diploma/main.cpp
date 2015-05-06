//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#include <stdio.h>

#include "field.h"

int main(int argc, const char * argv[]) {
    Field field;
    field.test();
    
    field.fillInitial();

    while (field.done() == false) {
        field.solve();
    }

    return 0;
}

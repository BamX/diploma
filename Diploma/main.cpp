//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#include <stdio.h>

#include "field.h"

int main(int argc, const char * argv[]) {
    Field field(6, 6, 100);
    field.fillInitial();

    while (field.done() == false) {
        field.solve();
        field.print();
    }

    return 0;
}

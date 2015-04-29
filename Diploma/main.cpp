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

    while (field.done() == false) {
        field.solve();

        field.print();

        if (field.time() > 5.0) {
            break;
        }
    }

    return 0;
}

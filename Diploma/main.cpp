//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#include <stdio.h>

#include "field.h"

int main(int argc, const char * argv[]) {
    Field field(100, 100, 20000, 0.01);
    field.test();

    field.enablePlotOutput();
    field.enableMatrixOutput();
    
    field.fillInitial();
    field.print();

    while (field.done() == false) {
        field.solve();

        field.print();

        if (field.time() > 1200) {
            break;
        }
    }

    return 0;
}

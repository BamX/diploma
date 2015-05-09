//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#include <stdio.h>
#include <mpi.h>

#include "field.h"

int main(int argc, char * argv[]) {
    MPI_Init(&argc, &argv);

    Field field;
    
    field.fillInitial();

    while (field.done() == false) {
        field.solve();

        if (field.time() > 5) {
            break;
        }
    }

    MPI_Finalize();
    return 0;
}

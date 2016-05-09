//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#include <stdio.h>
#include <mpi.h>

#include "field-static.h"
#include "field-transpose.h"
#include "factors.h"

int main(int argc, char * argv[]) {
    MPI_Init(&argc, &argv);

    {
        FieldTranspose field;
        //FieldStatic field;
        
        field.init();

        {
            const auto startTime = std::clock();

            while (field.done() == false) {
                field.solve();
            }

            const auto endTime = std::clock();

            int rank = 0;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            if (rank == 0) {
                std::cerr << "time: " << double(endTime - startTime) / CLOCKS_PER_SEC << '\n';
            }
        }

        field.finalize();
    }
    
    MPI_Finalize();
    return 0;
}

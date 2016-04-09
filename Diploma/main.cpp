//
//  Copyright (c) 2015 Nikolay Volosatov. All rights reserved.
//

#include <stdio.h>
#include <mpi.h>

#include "field.h"
#include "balancing.h"
#include "factors.h"

int ___main(int argc, char * argv[]) {
    std::vector<double> ww = {
        3.0, 40.0, 17.0, //  60
        10.0, 5.0, 17.0, 11.0, 11.0, 3.0, 1.0, // 58
        15.0, 17.0, 23.0, 6.0, 5.0, // 66
        66, // 66
        1.0, 1.0, 1.0, 1.0, 1.0, 62 // 67
    };

    size_t count = 5;
    auto parts = balancing::partition(ww, count);

    size_t i = 0;
    for (auto &part : parts) {
        double sum = 0;
        for (int j = 0; j < part; ++j) {
            std::cerr << ww[j + i] << " ";
            sum += ww[j + i];
        }
        i += part;
        std::cerr << "= " << sum << std::endl;
    }
    return 0;
}

int main(int argc, char * argv[]) {
    MPI_Init(&argc, &argv);

    Field field;
    
    field.fillInitial();

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
    MPI_Finalize();
    return 0;
}

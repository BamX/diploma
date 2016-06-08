//
//  Copyright © 2016 Nikolay Volosatov. All rights reserved.
//

#ifndef __Diploma__balancing__
#define __Diploma__balancing__

#include <vector>
#include <cstring>
#include <cmath>
#include <iostream>

namespace balancing {

    /**
     *  Weight-base partition.
     *
     *  @param weights job weights for partitioning
     *  @param lengths result intervals length
     */
    std::vector<int> partition(double *weights, size_t size, size_t count);

    std::vector<int> partition(std::vector<double> &weights, size_t count);

    std::vector<int> fastPartition(double *weights, size_t size, size_t *current_parts, size_t count);

    std::vector<int> fastPartition(std::vector<double> &weights, size_t *current_parts, size_t count);
    
}

#endif /* defined(__Diploma__balancing__)  */

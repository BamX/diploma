//
//  Copyright © 2016 Nikolay Volosatov. All rights reserved.
//

#ifndef __Diploma__balancing__
#define __Diploma__balancing__

#include <vector>

namespace balancing {

    std::vector<double> sums;

    void formatSums(double *weights, size_t size) {
        sums.resize(size);
        double prev = sums[0] = weights[0];
        for (size_t i = 1; i < size; ++i) {
            prev = sums[i] = prev + weights[i];
        }
    }

    size_t binSearch(size_t l, size_t r, double S) {
        size_t start = l;
        if (start > 0) {
            S += sums[start - 1];
        }

        while (true) {
            size_t mid = l + (r - l) / 2;

            if (sums[mid] > S) {
                if (l >= mid) {
                    return start;
                }

                r = mid;
            } else {
                if (mid == r - 1 || sums[mid + 1] > S) {
                    return mid + 1;
                }

                l = mid;
            }
        }
    }

    ssize_t bucketsCount(size_t start, size_t end, double S) {
        ssize_t partsCount = 0;
        while (start < end) {
            size_t pos = binSearch(start, end, S);
            if (pos == start) {
                return -1;
            }
            start = pos;
            partsCount += 1;
        }
        return partsCount;
    }

    bool buckets(size_t start, size_t end, double S, std::vector<size_t> &parts) {
        while (start < end) {
            size_t pos = binSearch(start, end, S);
            //if (start == pos) {
            //    return false;
            //}
            parts.push_back(pos - start);
            start = pos;
        }
        return true;
    }

    /**
     *  Weight-base partition.
     *
     *  @param weights job weights for partitioning
     *  @param lengths result intervals length
     */
    std::vector<size_t> partition(double *weights, size_t size, size_t count) {
        int i = 0, j = 0, iB = 0, jB = 0;
        double SB = __DBL_MAX__;

        formatSums(weights, size);

        while (i < size && j < size) {
            double S = sums[j] - (i > 0 ? sums[i - 1] : 0);
            ssize_t leftParts = bucketsCount(0, i, S);
            ssize_t rigthParts = bucketsCount(j + 1, size, S);
            if (leftParts < 0 || rigthParts < 0 || leftParts + rigthParts > count - 1) {
                ++j;
            } else {
                if (S < SB) {
                    iB = i; jB = j; SB = S;
                }
                if (i == j) {
                    ++j;
                }
                ++i;
            }
        }

        std::vector<size_t> parts;
        if (SB < __DBL_MAX__ - __DBL_EPSILON__) {
            parts.reserve(count);
            buckets(0, iB, SB, parts);
            parts.push_back(jB - iB + 1);
            buckets(jB + 1, size, SB, parts);
        }

        return parts;
    }

    std::vector<size_t> partition(std::vector<double> &weights, size_t count) {
        size_t size = weights.size();
        double *weightsAr = &weights[0];
        return partition(weightsAr, size, count);
    }
    
}

void test() {
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
}

#endif /* defined(__Diploma__balancing__)  */

//
//  Copyright © 2016 Nikolay Volosatov. All rights reserved.
//

#ifndef __Diploma__balancing__
#define __Diploma__balancing__

#include <vector>

namespace balancing {

    std::vector<double> sums;

    void formatSums(const std::vector<double> &weights) {
        sums.resize(weights.size());
        sums[0] = weights[0];
        for (size_t i = 1, len = weights.size(); i < len; ++i) {
            sums[i] = sums[i - 1] + weights[i];
        }
    }

    size_t binSearch(size_t l, size_t r, double S) {
        ssize_t start = l;
        double SStart = start > 0 ? sums[start - 1] : 0;

        while (true) {
            size_t mid = l + (r - l) / 2;
            double SMid = sums[mid] - SStart;

            if (SMid > S) {
                if (l + 1 >= mid) {
                    return start;
                }

                r = mid - 1;
            } else {
                if (mid == r - 1 || sums[mid + 1] - SStart > S) {
                    return mid + 1;
                }

                l = mid;
            }
        }
    }

    ssize_t bucketsCount(ssize_t start, ssize_t end, double S) {
        if (start >= end) {
            return 0;
        }

        ssize_t partsCount = 0;
        while (start < end) {
            ssize_t pos = binSearch(start, end, S);
            if (start == pos) {
                return -1;
            }
            start = pos;
            partsCount += 1;
        }
        return partsCount;
    }

    bool buckets(ssize_t start, ssize_t end, double S, std::vector<size_t> &parts) {
        if (start >= end) {
            return true;
        }

        while (start < end) {
            ssize_t pos = binSearch(start, end, S);
            if (start == pos) {
                return false;
            }
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
    std::vector<size_t> partition(const std::vector<double> &weights, size_t count) {
        size_t size = weights.size();

        int i = 0, j = 0, iB = 0, jB = 0;
        double SB = __DBL_MAX__;

        formatSums(weights);

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
        bool ok = buckets(0, iB, SB, parts);
        if (ok) {
            parts.push_back(jB - iB + 1);
        }
        ok = ok && buckets(jB + 1, size, SB, parts);

        return parts;
    }


    
}

#endif /* defined(__Diploma__balancing__)  */

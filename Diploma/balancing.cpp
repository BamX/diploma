//
//  balancing.cpp
//  Diploma
//
//  Created by Nikolay Volosatov on 08.06.16.
//  Copyright © 2016 Nikolay Volosatov. All rights reserved.
//

#include "balancing.h"

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

    long bucketsCount(size_t start, size_t end, double S) {
        long partsCount = 0;
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

    bool buckets(size_t start, size_t end, double S, std::vector<int> &parts) {
        while (start < end) {
            size_t pos = binSearch(start, end, S);
            //if (start == pos) {
            //    return false;
            //}
            parts.push_back((int)(pos - start));
            start = pos;
        }
        return true;
    }

    void fillNA(std::vector<int> &parts) {
        bool fixed = false;
        while (fixed == false) {
            fixed = true;
            bool better = false;
            for (size_t i = 0; i < parts.size(); ++i) {
                if (parts[i] == 0) {
                    if (i <= 1) {
                        if (parts[i + 1] > 1) {
                            parts[i + 1] -= 1;
                            parts[i] += 1;
                            better = true;
                        }
                    } else {
                        if (parts[i - 1] > 1) {
                            parts[i - 1] -= 1;
                            parts[i] += 1;
                            better = true;
                        }
                    }
                    fixed = false;
                }
            }
            if (fixed || (fixed == false && better == false)) {
                return;
            }
        }
    }

    /**
     *  Weight-base partition.
     *
     *  @param weights job weights for partitioning
     *  @param lengths result intervals length
     */
    std::vector<int> partition(double *weights, size_t size, size_t count) {
        int i = 0, j = 0, iB = 0, jB = 0;
        double SB = __DBL_MAX__;

        formatSums(weights, size);

        while (i < size && j < size) {
            double S = sums[j] - (i > 0 ? sums[i - 1] : 0);
            long leftParts = bucketsCount(0, i, S);
            long rigthParts = bucketsCount(j + 1, size, S);
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

        std::vector<int> parts;
        if (SB < __DBL_MAX__ - __DBL_EPSILON__) {
            parts.reserve(count);
            buckets(0, iB, SB, parts);
            parts.push_back(jB - iB + 1);
            buckets(jB + 1, size, SB, parts);
        }

        while (parts.size() < count) {
            parts.push_back(0);
        }
        fillNA(parts);

        return parts;
    }

    std::vector<int> partition(std::vector<double> &weights, size_t count) {
        size_t size = weights.size();
        double *weightsAr = &weights[0];
        return partition(weightsAr, size, count);
    }

    double costParts(std::vector<int> parts) {
        double SB = 0;
        double beginSum = 0;
        size_t beginIndex = 0;
        for (size_t i = 0; i < parts.size(); ++i) {
            double weightParts = 0;
            if (parts[i] > 0) {
                weightParts = sums[beginIndex + parts[i] - 1] - beginSum;
                beginSum += weightParts;
                beginIndex += parts[i];
            }
            if (weightParts > SB) {
                SB = weightParts;
            }
        }
        return SB;
    }

    std::vector<int> fastPartition(double *weights, size_t size, size_t *current_parts, size_t count) {
        double SB = __DBL_MAX__;

        formatSums(weights, size);

        std::vector<int> parts(current_parts, current_parts + count);

        size_t iterations = 0;
        size_t iterationsLimit = log(size);
        while (iterations < iterationsLimit) {
            double newSB = costParts(parts);
            if (fabs(SB - newSB) < 0.00001) {
                break;
            }
            ++iterations;
            SB = newSB;
            double S = 0;
            size_t idx = 0;
            for (size_t i = 0; i < count - 1; ++i) {
                double s1 = sums[idx + parts[i] - 1] - S;
                double s2 = sums[idx + parts[i] + parts[i + 1] - 1] - sums[idx + parts[i] - 1];
                double z1 = weights[idx + parts[i] - 1];
                double z2 = weights[idx + parts[i]];
                if (s1 + z2 < SB && s1 < s2) {
                    parts[i] += 1;
                    parts[i + 1] -= 1;
                } else if (s2 + z1 < SB && s1 > s2) {
                    parts[i + 1] += 1;
                    parts[i] -= 1;
                }

                S = sums[idx + parts[i] - 1];
                idx += parts[i];
            }
        }

        return parts;
    }

    std::vector<int> fastPartition(std::vector<double> &weights, size_t *current_parts, size_t count) {
        size_t size = weights.size();
        double *weightsAr = &weights[0];
        return fastPartition(weightsAr, size, current_parts, count);
    }

}

int test_main() {
    std::vector<double> ww = {
        3.0, 40.0, 17.0, //  60
        10.0, 5.0, 17.0, 11.0, 11.0, 3.0, 1.0, // 58
        15.0, 17.0, 23.0, 6.0, 5.0, // 66
        66, // 66
        1.0, 1.0, 1.0, 1.0, 1.0, 62 // 67
    };
    size_t lastParts[] = { 4, 4, 5, 5, 4 };

    size_t count = 5;
    auto parts = balancing::fastPartition(ww, lastParts, count);

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

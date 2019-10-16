// Author: Thomas Miethlinger BSc.
// thomas.miethlinger@gmail.com

// STL
#include <vector>
#include <algorithm>
#include <functional>

// Boost
#include <boost/algorithm/string.hpp>

#pragma once

typedef std::vector<int> VI;
typedef std::vector<VI> VVI;
typedef std::vector<double> VD;
typedef std::vector<VD> VVD;
using std::size_t;

namespace RQA
{
    VD CalculateRecurrenceMatrix(VD &norm_distance_matrix, double eps)
    {
        VD recurrence_matrix(norm_distance_matrix.size(), 0);
        std::transform(norm_distance_matrix.begin(), norm_distance_matrix.end(), recurrence_matrix.begin(), [eps](double d){ return (d <= eps ? 1 : 0); };
        return recurrence_matrix;
    }

    VD NormalizeDistanceMatrix(VD &distance_matrix)
    {
        VD norm_distance_matrix(distance_matrix.size(), 0);
        double min = *std::min_element(distance_matrix.begin(), distance_matrix.end());
        double max = *std::max_element(distance_matrix.begin(), distance_matrix.end());
        std::transform(distance_matrix.begin(), distance_matrix.end(), norm_distance_matrix.begin(), [min, max](double d){ return ((d - min)/(max - min)); };
        return norm_distance_matrix;
    }
};

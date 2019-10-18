// Author:     Thomas Miethlinger B.Sc.
// Date:       16.10.2019
// Email:      t.miethlinger@protonmail.com
// University: Johannes Kepler University Linz, Austria

// STL
#include<vector>
#include<list>
#include<cmath> // sqrt
#include<algorithm> // std::transform, std::min_element, std::max_element, std::fill
#include<numeric> // std::accumulate

#pragma once

namespace rqa_core
{
    void calculate_recurrence_matrix(std::vector<double> &norm_distance_matrix, double eps, std::vector<int> &recurrence_matrix)
    {
        std::transform(norm_distance_matrix.begin(), norm_distance_matrix.end(), recurrence_matrix.begin(), [eps](double d){ return (d <= eps ? 1 : 0); });
    }

    void normalize_distance_matrix(std::vector<double> &distance_matrix, std::vector<double> &norm_distance_matrix)
    {
        double min = *std::min_element(distance_matrix.begin(), distance_matrix.end());
        double max = *std::max_element(distance_matrix.begin(), distance_matrix.end());
        std::transform(distance_matrix.begin(), distance_matrix.end(), norm_distance_matrix.begin(), [min, max](double d){ return ((d - min)/(max - min)); });
    }

    // Functions for recurrence quantitative analysis (RQA)
    // Based on http://www.recurrence-plot.tk/rqa.php

    // Calculate the recurrence rate rr for a given recurrence matrix
    double recurrence_rate(std::vector<int> &recurrence_matrix)
    {
        double rr = 0;
        int n2 = recurrence_matrix.size();
        double sum = std::accumulate(recurrence_matrix.begin(), recurrence_matrix.end(), 0.0);
        rr = sum / (double)n2;
        return rr;
    }

    // Calculate the determinism vector for a given recurrence matrix
    std::vector<double> determinism(std::vector<int> &recurrence_matrix)
    {
        // determinism vector is of size n due to possible min and max
        // sizes of diagonal lines l
        int n2 = (int)recurrence_matrix.size();
        int n = sqrt(n2);
        std::vector<double> det_vector(n, 0);

        // Create vector which represents histogram of lengths of diagonal lines of size n:
        // Smallest possible diagonal line size = 1
        // Largest possible diagonal line size = n -> vector.size == n
        std::vector<int> hist_diagonal_lines(n, 0);
        // There is always a recurrence in the main diagonal
        hist_diagonal_lines[n - 1] = 1;

        // Iterate through the half of the recurrence_matrix in an diagonal manner (without main diagonal)
        std::size_t row;
        std::size_t col;
        bool found_diag;
        int ctr_diag_size;
        for(int i = 0; i < n; i++)
        {
            found_diag = false;
            for(int j = 0; j < i + 1; j++)
            {
                row = n - i + j - 1;
                col = j;

                // Find diagonal lines and save them in the vector hist_diagonal_lines
                if(recurrence_matrix[row * n + col] == 1 && !found_diag)
                {
                    found_diag = true;
                    ctr_diag_size = 1;
                }
                else if(recurrence_matrix[row * n + col] == 1 && found_diag)
                {
                    ctr_diag_size++;
                }
                else if(recurrence_matrix[row * n + col] == 0 && found_diag)
                {
                    found_diag = false;
                    // There are always two diagonal lines; top and bottom of diagonal lines
                    hist_diagonal_lines[ctr_diag_size - 1] += 2;
                }
                // else: recurrence_matrix[row * n + col] == 0 && found_diag -> nothing to do
            }
        }

        // Compute the determinism vector

        // Denominator is independent of lmin
        double denom = 0.0;
        for(int l = 1; l <= n; l++)
        {
            // Histogram of diagonal line of length l is stored
            // in hist_diagonal_lines at index l - 1
            denom += l * hist_diagonal_lines[l - 1];
        }

        // Numerator is dependent of lmin
        double num;
        for(int lmin = 1; lmin <= n; lmin++)
        {
            num = 0.0;
            for(int l = lmin; l <= n; l++)
            {
                num += l * hist_diagonal_lines[l - 1];
            }
            det_vector[lmin - 1] = num / denom;
        }

        return det_vector;
    }

    // Calculate the laminarity vector for a given recurrence matrix
    std::vector<double> laminarity(std::vector<int> &recurrence_matrix)
    {
        // laminarity vector is of size n due to possible min and max
        // sizes of vertical lines l
        int n2 = (int)recurrence_matrix.size();
        int n = sqrt(n2);
        std::vector<double> lam_vector(n, 0);

        // Create vector which represents histogram of lengths of vertical(and, by mirroring, horizontal) lines of size n:
        // Smallest possible vertical line size = 1
        // Largest possible vertical line size = n -> vector.size == n
        std::vector<int> hist_vertical_lines(n, 0);

        // Iterate through the half of the recurrence_matrix in an vertical manner (without main vertical)
        std::size_t row;
        std::size_t col;
        bool found_vert;
        int ctr_vert_size;
        for(int i = 0; i < n; i++)
        {
            found_vert = false;
            for(int j = 0; j < i + 1; j++)
            {
                row = i - j;
                col = i;

                // Find vertical lines and save them in the vector hist_vertical_lines
                if(recurrence_matrix[row * n + col] == 1 && !found_vert)
                {
                    found_vert = true;
                    ctr_vert_size = 1;
                }
                else if(recurrence_matrix[row * n + col] == 1 && found_vert)
                {
                    ctr_vert_size++;
                }
                else if(recurrence_matrix[row * n + col] == 0 && found_vert)
                {
                    found_vert = false;
                    // There are always two vertical lines; top and bottom of vertical lines
                    hist_vertical_lines[ctr_vert_size - 1] += 2;
                }
                // else: recurrence_matrix[row * n + col] == 0 && found_vert -> nothing to do
            }
        }

        // Compute the laminarity vector

        // Denominator is independent of lmin
        double denom = 0.0;
        for(int l = 1; l <= n; l++)
        {
            // Histogram of vertical line of length l is stored
            // in hist_vertical_lines at index l - 1
            denom += l * hist_vertical_lines[l - 1];
        }

        // Numerator is dependent of lmin
        double num;
        for(int lmin = 1; lmin <= n; lmin++)
        {
            num = 0.0;
            for(int l = lmin; l <= n; l++)
            {
                num += l * hist_vertical_lines[l - 1];
            }
            lam_vector[lmin - 1] = num / denom;
        }

        return lam_vector;
    }
};

// Author:     Thomas Miethlinger B.Sc.
// Date:       21.10.2019
// Email:      t.miethlinger@protonmail.com
// University: Johannes Kepler University Linz, Austria

#include <algorithm> // std::nth_element
#include <functional> // std::function
#include <utility> // std::pair
#include <vector> // std::vector
#include <list> // std::list

#pragma once

namespace ha_cost
{
    using std::function;
    using std::pair;
    using std::vector;

    typedef vector<int> VI;
    typedef vector<double> VD;
    typedef vector<VD> VVD;
    typedef pair<int, double> PID;

    // Maximum cost, symbolizing infinity
    constexpr double cost_max = (double)std::numeric_limits<int>::max();

    VVD create_costobject_adjmatrix(int N, VVD &Xa, VVD &Xb, function<double(VD&, VD&, double, double)> d, double l_inv, double v_inv)
    {
        // Define cost matrix object
        VVD cost_adjmatrix(N, VD(N, 0.0));

        // Compute cost matrix for each pair (i, j)
        // from the distance d
        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < N; j++)
            {
                cost_adjmatrix[i][j] = d(Xa[i], Xb[j], l_inv, v_inv);
            }
        }

        return cost_adjmatrix;
    }

    VVD create_costobject_adjmatrix(int N, int N_neighbours, VVD &Xa, VVD &Xb, function<double(VD&, VD&, double, double)> d, double l_inv, double v_inv)
    {
        // Define cost matrix object with standard value cost_max
        VVD cost_adjmatrix(N, VD(N, cost_max));

        // Compute all distances and set the N_neighbours smallest ones for each row
        vector<PID> cost_node_row(N);
        // For each row i do:
        for(int i = 0; i < N; i++)
        {
            // For each column j do:
            // Compute the distance d_ij and store it together with the
            // corresponding index j into the j-th entry of cost_node_row
            for(int j = 0; j < N; j++)
            {
                PID p;
                p.first = j;
                p.second = d(Xa[i], Xb[j], l_inv, v_inv);
                cost_node_row[j] = p;
            }

            // Partial sort up to the n-th element of cost_node_row
            std::nth_element(
                                cost_node_row.begin(), 
                                cost_node_row.begin() + N_neighbours, 
                                cost_node_row.end(), 
                                [](PID a, PID b) {return a.second < b.second;}
                            );

            // Assign N_neighbours smallest values in cost_adjmatrix
            for(int j = 0; j < N_neighbours; j++)
                cost_adjmatrix[i][cost_node_row[j].first] = cost_node_row[j].second;
        }

        // Do the same column-wise
        vector<PID> cost_node_col(N);
        for(int j = 0; j < N; j++)
        {
            for(int i = 0; i < N; i++)
            {
                PID p;
                p.first = i;
                p.second = d(Xa[i], Xb[j], l_inv, v_inv);
                cost_node_col[i] = p;
            }

            std::nth_element(
                                cost_node_col.begin(), 
                                cost_node_col.begin() + N_neighbours, 
                                cost_node_col.end(), 
                                [](PID a, PID b) {return a.second < b.second;}
                            );

            for(int i = 0; i < N_neighbours; i++)
                cost_adjmatrix[cost_node_col[i].first][j] = cost_node_col[i].second;
        }

        return cost_adjmatrix;
    }

    vector<vector<PID>> create_costobject_adjlist_plain(int N, int N_neighbours, VVD &Xa, VVD &Xb, function<double(VD&, VD&, double, double)> d, double l_inv, double v_inv)
    {
        // Calculate
        VVD cost_adjmatrix = create_costobject_adjmatrix(N, N_neighbours, Xa, Xb, d, l_inv, v_inv);

        // Define cost matrix object with standard value cost_max
        vector<vector<PID>> cost_adjlist(N);

        int ctr;
        for(int i = 0; i < N; i++)
        {
            ctr = 0;
            for(int j = 0; j < N; j++)
            {
                ctr = cost_adjmatrix[i][j] < cost_max ? (ctr + 1) : ctr;
            }
            cost_adjlist[i].resize(ctr);

            for(int j = 0, k = 0; j < N; j++)
            {
                if(cost_adjmatrix[i][j] < cost_max)
                {
                    cost_adjlist[i][k] = std::make_pair(j, cost_adjmatrix[i][j]);
                    k++;
                }
            }
        }

        return cost_adjlist;
    }

    vector<vector<PID>> create_costobject_adjlist_slim(int N, int N_neighbours, VVD &Xa, VVD &Xb, function<double(VD&, VD&, double, double)> d, double l_inv, double v_inv)
    {
        // Define cost matrix object with standard value cost_max
        vector<vector<PID>> cost_adjlist(N, vector<PID>(2 * N_neighbours));

        // Find the M nearest neighbours per row
        vector<vector<PID>> cost_adjlist_row(N, vector<PID>(N_neighbours));
        vector<PID> cost_node_row(N);
        vector<PID> cost_node_row_part(N_neighbours);
        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < N; j++)
            {
                cost_node_row[j] = std::make_pair(j, d(Xa[i], Xb[j], l_inv, v_inv));
            }

            std::nth_element(cost_node_row.begin(), cost_node_row.begin() + N_neighbours, cost_node_row.end(), [](PID a, PID b) {return a.second < b.second;});
            std::copy(cost_node_row.begin(), cost_node_row.begin() + N_neighbours, cost_node_row_part.begin());
            std::sort(cost_node_row_part.begin(), cost_node_row_part.end(), [](PID a, PID b) {return a.second < b.second;});
            std::copy(cost_node_row_part.begin(), cost_node_row_part.end(), cost_adjlist_row[i].begin());
        }

        // Find the M nearest neighbours per column
        vector<vector<PID>> cost_adjlist_col(N, vector<PID>(N_neighbours));
        vector<PID> cost_node_col(N);
        vector<PID> cost_node_col_part(N_neighbours);
        for(int j = 0; j < N; j++)
        {
            for(int i = 0; i < N; i++)
            {
                cost_node_col[i] = std::make_pair(i, d(Xa[i], Xb[j], l_inv, v_inv));
            }

            std::nth_element(cost_node_col.begin(), cost_node_col.begin() + N_neighbours, cost_node_col.end(), [](PID a, PID b) {return a.second < b.second;});
            std::copy(cost_node_col.begin(), cost_node_col.begin() + N_neighbours, cost_node_col_part.begin());
            std::sort(cost_node_col_part.begin(), cost_node_col_part.end(), [](PID a, PID b) {return a.second < b.second;});
            std::copy(cost_node_col_part.begin(), cost_node_col_part.end(), cost_adjlist_col[j].begin());
        }

        // Merge such that the number of neighbours is at least m row-&columnwise
        vector<std::list<PID>> cost_adjlist_tmp(N, std::list<PID>());
        for(int i = 0; i < N; i++)
        {
            cost_adjlist_tmp[i].assign(cost_adjlist_row[i].begin(), cost_adjlist_row[i].end());
        }
        bool found;
        for(int j = 0; j < N; j++)
        {
            for(vector<PID>::iterator i_it = cost_adjlist_col[j].begin(); i_it != cost_adjlist_col[j].end(); i_it++)
            {
                int i = i_it->first;
                found = false;
                for(vector<PID>::iterator j_it = cost_adjlist_row[i].begin(); j_it != cost_adjlist_row[i].end(); j_it++)
                {
                    if(j == j_it->first)
                    {
                        found = true;
                        break;
                    }
                }
                if(!found)
                {
                    cost_adjlist_tmp[i].push_back(std::make_pair(j, i_it->second));
                }
            }
        }

        // Sort the merged lists
        for(int i = 0; i < N; i++)
        {
            cost_adjlist_tmp[i].sort([](PID a, PID b) {return a.first < b.first;});
        }

        // Assign lists from cost_adjlist_tmp to vectors of cost_adjlist
        for(int i = 0; i < N; i++)
        {
            cost_adjlist[i].assign(cost_adjlist_tmp[i].begin(), cost_adjlist_tmp[i].end());
        }

        return cost_adjlist;
    }

    // Compute costs from given matching solution
    double compute_total_cost_adjmatrix(int N, VVD& cost_adjmatrix, VI q)
    {
        double total_cost = 0.0;
        int count = 0;

        for(int i = 0; i < N; i++)
        {
            if(q[i] != -1)
            {
                count++;
                total_cost += cost_adjmatrix[i][q[i]];
            }
        }

        if(count != N)
        {
            total_cost *= ((double) N / (double) count);
        }

        return total_cost;
    }

    double compute_total_cost_adjlist(int N, vector<vector<PID>>& cost_adjlist, VI q)
    {
        double total_cost = 0.0;
        int count = 0;
        int size;

        for(int i = 0; i < N; i++)
        {
            if(q[i] != -1)
            {
                size = cost_adjlist[i].size();
                bool found = false;
                for(int j = 0; j < size; j++)
                {
                    if(q[i] == cost_adjlist[i][j].first)
                    {
                        total_cost += cost_adjlist[i][j].second;
                        count++;
                        found = true;
                    }
                }
                if(!found)
                    total_cost += cost_max;
            }
        }

        if(count != N)
        {
            total_cost *= ((double) N / (double) count);
        }

        return total_cost;
    }
};

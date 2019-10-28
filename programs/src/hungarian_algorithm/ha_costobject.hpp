// Author:     Thomas Miethlinger B.Sc.
// Date:       21.10.2019
// Email:      t.miethlinger@protonmail.com
// University: Johannes Kepler University Linz, Austria

#include <algorithm> // std::nth_element
#include <functional> // std::function
#include <utility> // std::pair
#include <vector> // std::vector

#pragma once

namespace ha_costobject
{
    using std::function;
    using std::pair;
    using std::vector;

    typedef vector<double> VD;
    typedef vector<VD> VVD;
    typedef pair<int, double> PID;

    // Maximum cost, symbolizing infinity
    constexpr double cost_max = (double)std::numeric_limits<int>::max();

    VVD create_costobject_adjmatrix(int N, VVD &Xa, VVD &Xb, function<double(VD&, VD&, double, double)> d, double l, double v)
    {
        // Define cost matrix object
        VVD cost_adjmatrix(N, VD(N, 0.0));

        // Invert l and v and store result to obtain speed-up
        double l_inv = 1.0 / l;
        double v_inv = 1.0 / v;

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

    VVD create_costobject_adjmatrix(int N, int N_neighbours, VVD &Xa, VVD &Xb, function<double(VD&, VD&, double, double)> d, double l, double v)
    {
        // Define cost matrix object with standard value cost_max
        VVD cost_adjmatrix(N, VD(N, cost_max));

        // Invert l and v and store result to obtain speed-up
        double l_inv = 1.0 / l;
        double v_inv = 1.0 / v;

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

    vector<vector<PID>> create_costobject_adjlist_plain(int N, int N_neighbours, VVD &Xa, VVD &Xb, function<double(VD&, VD&, double, double)> d, double l, double v)
    {
        // Calculate
        VVD cost_adjmatrix = create_costobject_adjmatrix(N, N_neighbours, Xa, Xb, d, l, v);

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

    vector<vector<PID>> create_costobject_adjlist_slim(int N, int N_neighbours, VVD &Xa, VVD &Xb, function<double(VD&, VD&, double, double)> d, double l, double v)
    {
        // Define cost matrix object with standard value cost_max
        vector<vector<PID>> cost_adjlist(N, vector<PID>(2 * N_neighbours));

        // Invert l and v and store result to obtain speed-up
        double l_inv = 1.0 / l;
        double v_inv = 1.0 / v;

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
};

    // Create cost object for assignment using hungarian algorithm
    // create_co ... create cost object
    // amatrix ... use full adjacency matrix as representation for cost graph
    // alist ... use adjacency list as representation for cost graph
    // nearnb ... make cost graph sparse by only considering the M nearest neighbours

    // Finite Domain Boundary Conditions
    //public:

    // VVD create_cost_matrix_naive_position_nearnb(VVD &Ra, VVD &Rb, int M);
    // vector<vector<ParticleAssignment::SkeletonNode>> create_cost_matrix_adjlist_position_nearnb(VVD &Ri, VVD &Rj, int M);






    /*static vector<vector<double>> CreateCostMatrix(vector<vector<double>> &Ri, vector<vector<double>> &Vi,
    vector<vector<double>> &Rj, vector<vector<double>> &Vj);
    static vector<vector<double>> CreateCostMatrix(vector<vector<double>> &Ri, vector<vector<double>> &Vi,
    vector<vector<double>> &Rj, vector<vector<double>> &Vj, int M);
    static vector<vector<ParticleAssignment::SkeletonNode>> CreateCostSkeleton(vector<vector<double>> &Ri, vector<vector<double>> &Vi,
    vector<vector<double>> &Rj, vector<vector<double>> &Vj, int M);
    static vector<vector<double>> CreateCostMatrixPBC(vector<vector<double>> &Ri, vector<vector<double>> &Vi,
    vector<vector<double>> &Rj, vector<vector<double>> &Vj, double L, double L2);*/

/*

vector<vector<double>> ParticleAssignment::CreateCostMatrix(vector<vector<double>> &Ri, vector<vector<double>> &Rj)
{
    int N = Ri.size();
    vector<vector<double>> costmatrix(N, vector<double>(N));

    for(int i = 0; i < N; i++)
    {
    for(int j = 0; j < N; j++)
    {
        costmatrix[i][j] = EuclideanMetric(Ri[i], Rj[j]);
    }
    }

    return costmatrix;
}

class ParticleAssignment
{
    public:
    static ParticleAssignment::Result MinCostMatching(vector<vector<double>> &Ri, vector<vector<double>> &Rj, int M, int B);
    static ParticleAssignment::Result MinCostMatchingMatrix(vector<vector<double>> &costmatrix);
    static ParticleAssignment::Result MinCostMatchingSkeleton(vector<vector<ParticleAssignment::SkeletonNode>> &costskeleton);
    static vector<vector<double>> CreateCostMatrix(vector<vector<double>> &Ri, vector<vector<double>> &Rj);
    static vector<vector<double>> CreateCostMatrix(vector<vector<double>> &Ri, vector<vector<double>> &Rj, int M);
    static vector<vector<ParticleAssignment::SkeletonNode>> CreateCostSkeleton(vector<vector<double>> &Ri, vector<vector<double>> &Rj, int M);


    private:


    static ParticleAssignment::CostObject CostFromMatchingMatrix(vector<vector<double>> &costmatrix, vector<int> &q);
    static ParticleAssignment::CostObject CostFromMatchingSkeleton(vector<vector<ParticleAssignment::SkeletonNode>> &costskeleton, vector<int> &q);
    static ParticleAssignment::CostObject CostFromMatchingVectors(vector<vector<double>> &Ri, vector<vector<double>> &Rj, vector<int> &q);

    static double EuclideanMetric(vector<double> &ri, vector<double> &rj);

    static constexpr double precision = 0.000000001;
    static constexpr double cost_max = (double)std::numeric_limits<int>::max();
};*/

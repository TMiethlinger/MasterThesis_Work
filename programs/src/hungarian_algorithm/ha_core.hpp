// Author:     Thomas Miethlinger B.Sc.
// Date:       21.10.2019
// Email:      t.miethlinger@protonmail.com
// University: Johannes Kepler University Linz, Austria

#include <cmath>
#include <vector>
#include <limits>

#pragma once

namespace ha_core
{
    using std::vector;

    typedef vector<int> VI;
    typedef vector<double> VD;
    typedef vector<VD> VVD;

    // Precision which is used to check for zero cost
    constexpr double precision = 0.000000001;

    // Maximum cost, symbolizing infinity
    constexpr double cost_max = (double)std::numeric_limits<int>::max();


    // Compute the hungarian algorithm minimum cost matching solution
    // for an adjacency cost matrix
    VI min_cost_matching_adjmatrix(VVD &cost_adjmatrix)
    {
        int N = cost_adjmatrix.size();

        // construct dual feasible solution
        VD u(N, cost_max);
        VD v(N, cost_max);
        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < N; j++)
            {
                // Minimize cost_adjmatrix for each row i w.r.t. all columns j
                u[i] = std::min(u[i], cost_adjmatrix[i][j]);
            }
        }
        for(int j = 0; j < N; j++)
        {
            for(int i = 0; i < N; i++)
            {
                // Minimize cost_adjmatrix - u for each column j w.r.t. all rows i
                v[j] = std::min(v[j], cost_adjmatrix[i][j] - u[i]);
            }
        }

        // construct primal solution satisfying complementary slackness
        VI q = VI(N, -1);
        VI q_inv = VI(N, -1);
        int mated = 0;
        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < N; j++)
            {
                // if q_inv[j] has a matching we can continue
                if (q_inv[j] != -1)
                    continue;
                // if cost_adjmatrix[i][j] - u[i] - v[j] == 0
                // we match (i, j) and increase the counter mated
                if (fabs(cost_adjmatrix[i][j] - u[i] - v[j]) < precision)
                {
                    q[i] = j;
                    q_inv[j] = i;
                    mated++;
                    break;
                }
            }
        }

        VD dist(N);
        VI dad(N);
        VI seen(N);

        // repeat until primal solution is feasible
        while (mated < N)
        {
            // find an unmatched left node
            // store the corresponding index of this single node in s
            int s = 0;
            while (q[s] != -1)
                s++;
 
            // initialize Dijkstra
            std::fill(dad.begin(), dad.end(), -1);
            std::fill(seen.begin(), seen.end(), 0);
            for(int k = 0; k < N; k++)
            {
                dist[k] = cost_adjmatrix[s][k] - u[s] - v[k];
            }

            int j = 0;
            while (true)
            {
                // find closest
                j = -1;
                for(int k = 0; k < N; k++)
                {
                    if (!seen[k])
                    {
                        if (j == -1 || dist[k] < dist[j])
                            j = k;
                    }
                }
                seen[j] = 1;

                // termination condition
                if (q_inv[j] == -1)
                    break;

                // relax neighbours
                const int i = q_inv[j];
                for(int k = 0; k < N; k++)
                {
                    if (!seen[k])
                    {
                        const double new_dist = dist[j] + cost_adjmatrix[i][k] - u[i] - v[k];
                        if (new_dist < dist[k])
                        {
                            dist[k] = new_dist;
                            dad[k] = j;
                        }
                    }
                }
            }

            // update dual variables
            for(int k = 0; k < N; k++)
            {
                if (k != j && seen[k])
                {
                    const double delta = dist[k] - dist[j];
                    v[k] += delta;
                    u[q_inv[k]] -= delta;
                }
            }
            u[s] += dist[j];

            // augment along path
            while (dad[j] >= 0)
            {
                const int d = dad[j];
                q_inv[j] = q_inv[d];
                q[q_inv[j]] = j;
                j = d;
            }
            q_inv[j] = s;
            q[s] = j;

            mated++;
        }

        return q;
    }
};

// Author:     Thomas Miethlinger B.Sc.
// Date:       21.10.2019
// Email:      t.miethlinger@protonmail.com
// University: Johannes Kepler University Linz, Austria

#include <cmath>
#include <vector>
#include <limits>
#include <utility> // std::pair

#pragma once

namespace ha_core
{
    using std::function;
    using std::list;
    using std::pair;
    using std::set;
    using std::vector;

    typedef vector<int> VI;
    typedef vector<double> VD;
    typedef vector<VD> VVD;
    typedef pair<int, double> PID;

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
                if(q_inv[j] != -1)
                    continue;
                // if cost_adjmatrix[i][j] - u[i] - v[j] == 0
                // we match (i, j) and increase the counter mated
                if(fabs(cost_adjmatrix[i][j] - u[i] - v[j]) < precision)
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
        while(mated < N)
        {
            // find an unmatched left node
            // store the corresponding index of this single node in s
            int s = 0;
            while(q[s] != -1)
                s++;
 
            // initialize Dijkstra
            std::fill(dad.begin(), dad.end(), -1);
            std::fill(seen.begin(), seen.end(), 0);
            for(int k = 0; k < N; k++)
            {
                dist[k] = cost_adjmatrix[s][k] - u[s] - v[k];
            }

            int j = 0;
            while(true)
            {
                // find closest
                j = -1;
                for(int k = 0; k < N; k++)
                {
                    if(!seen[k])
                    {
                        if(j == -1 || dist[k] < dist[j])
                            j = k;
                    }
                }
                seen[j] = 1;

                // termination condition
                if(q_inv[j] == -1)
                    break;

                // relax neighbours
                const int i = q_inv[j];
                for(int k = 0; k < N; k++)
                {
                    if(!seen[k])
                    {
                        const double new_dist = dist[j] + cost_adjmatrix[i][k] - u[i] - v[k];
                        if(new_dist < dist[k])
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
                if(k != j && seen[k])
                {
                    const double delta = dist[k] - dist[j];
                    v[k] += delta;
                    u[q_inv[k]] -= delta;
                }
            }
            u[s] += dist[j];

            // augment along path
            while(dad[j] >= 0)
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

    VI min_cost_matching_adjlist(vector<vector<PID>> &cost_adjlist)
    {
        int N = cost_adjlist.size();

        function<bool(PID, PID)> comparator_PID = 
        [](PID x1, PID x2)
        {
            if(x1.second < x2.second)
                return true;
            else if(x1.second > x2.second)
                return false;
            else
                return x1.first < x2.first;
         };

        // Construct dual feasible solution
        VD u(N, cost_max);
        VD v(N, cost_max);
        for(int i = 0; i < N; i++)
        {
            const int size = cost_adjlist[i].size();
            for(int j = 0; j < size; j++)
            {
                u[i] = std::min(u[i], cost_adjlist[i][j].second);
            }
        }
        for(int i = 0; i < N; i++)
        {
            const int size = cost_adjlist[i].size();
            for(int j = 0, y; j < size; j++)
            {
                y = cost_adjlist[i][j].first;
                v[y] = std::min(v[y], cost_adjlist[i][j].second - u[i]);
            }
        }

        // Construct primal solution satisfying complementary slackness
        VI q = VI(N, -1);
        VI q_inv = VI(N, -1);
        int mated = 0;
        for(int i = 0; i < N; i++)
        {
            const int size = cost_adjlist[i].size();
            for(int j = 0, y; j < size; j++)
            {
                y = cost_adjlist[i][j].first;
                if(q_inv[y] == -1)
                {
                    if(fabs(cost_adjlist[i][j].second - u[i] - v[y]) < precision)
                    {
                        q[i] = y;
                        q_inv[y] = i;
                        mated++;
                        break;
                    }
                }
            }
        }

        VD dist(N);
        VI dad(N);
        VI seen(N);

        // Memory variable for s: allows to skip non-matched nodes,
        // and go to the next one
        list<int> slist;
        for(int i = 0; i < N; i++)
            slist.push_back(i);

        // Repeat until primal solution is feasible
        int s = -1;
        list<int>::iterator iter_s;
        set<PID, function<bool(PID, PID)>> vertices_right_set(comparator_PID);
        pair<bool, set<PID, function<bool(PID, PID)>>::iterator> default_pair = std::make_pair(false, vertices_right_set.end());
        vector<pair<bool, set<PID, function<bool(PID, PID)>>::iterator>> active_vertices(N, default_pair);
        // repeat until primal solution is feasible
        while(mated < N)
        {
            // if: there was an assignment previously and consequently s == -1
            if(s == -1)
            {
                // Find an unmatched left node
                iter_s = slist.begin();
                while(q[*iter_s] != -1 && iter_s != slist.end())
                {
                    iter_s++;
                }
                if(iter_s == slist.end())
                {
                    // no other assignments free
                    return q;
                }
                else
                {
                    s = *iter_s;
                }
            }
            // else: Use previously found unmatched left node

            // Initialize Dijkstra
            std::fill(dad.begin(), dad.end(), -1);
            std::fill(seen.begin(), seen.end(), 0);
            std::fill(dist.begin(), dist.end(), cost_max);
            const int size = cost_adjlist[s].size();

            for(int l = 0, y; l < size; l++)
            {
                y = cost_adjlist[s][l].first;
                dist[y] = cost_adjlist[s][l].second - u[s] - v[y];
                if(active_vertices[y].first)
                {
                    // The index y is already active:
                    // 1. Delete existing entry
                    vertices_right_set.erase(active_vertices[y].second);
                    // 2. Add entry and set iterator
                    active_vertices[y].second = vertices_right_set.insert(std::make_pair(y, dist[y])).first;
                }
                else
                {
                    // The index y is inactive
                    // 1. Set index active
                    active_vertices[y].first = true;
                    // 2. Add entry and set iterator
                    active_vertices[y].second = vertices_right_set.insert(std::make_pair(y, dist[y])).first;
                }
            }

            int j;
            while(true)
            {
                j = -1;
                set<PID>::iterator iter_min = vertices_right_set.begin();
                if(iter_min != vertices_right_set.end())
                {
                    j = iter_min->first;
                    active_vertices[j].first = false;
                    active_vertices[j].second = vertices_right_set.end();
                    vertices_right_set.erase(iter_min);
                    seen[j] = 1;
                }
                else
                {
                    // No assignment for this s possible
                    // Assume that this will not be changed by later assignments
                    slist.erase(iter_s);

                    // search for new s
                    iter_s = slist.begin();
                    while(q[*iter_s] != -1 && iter_s != slist.end())
                    {
                        iter_s++;
                    }
                    if(iter_s == slist.end())
                    {
                        // no other assignments free
                        return q;
                    }
                    else
                    {
                        s = *iter_s;
                        break;
                    }
                }

                // Termination condition
                if(q_inv[j] == -1)
                    break;

                const int i = q_inv[j];
                const int size = cost_adjlist[i].size();
                for(int l = 0, y; l < size; l++)
                {
                    y = cost_adjlist[i][l].first;
                    if(!seen[y])
                    {
                        const double new_dist = dist[j] + cost_adjlist[i][l].second - u[i] - v[y];
                        if(new_dist < dist[y])
                        {
                            dist[y] = new_dist;
                            dad[y] = j;
                            if(active_vertices[y].first)
                            {
                                // The index y is already active:
                                // 1. Delete existing entry
                                vertices_right_set.erase(active_vertices[y].second);
                                // 2. Add entry and set iterator
                                active_vertices[y].second = vertices_right_set.insert(std::make_pair(y, dist[y])).first;
                            }
                            else
                            {
                                // The index y is inactive
                                // 1. Set index active
                                active_vertices[y].first = true;
                                // 2. Add entry and set iterator
                                active_vertices[y].second = vertices_right_set.insert(std::make_pair(y, dist[y])).first;
                            }
                        }
                    }
                }
            }

            if(j != -1)
            {
                // Update dual variables
                for(int k = 0; k < N; k++)
                {
                    if(k != j && seen[k])
                    {
                        const double delta = dist[k] - dist[j];
                        v[k] += delta;
                        u[q_inv[k]] -= delta;
                    }
                }
                u[s] += dist[j];

                // Augment along path
                while(dad[j] >= 0)
                {
                    const int d = dad[j];
                    q_inv[j] = q_inv[d];
                    q[q_inv[j]] = j;
                    j = d;
                }
                q_inv[j] = s;
                q[s] = j;

                mated++;
                vertices_right_set.clear();
                std::fill(active_vertices.begin(), active_vertices.end(), default_pair);
                s = -1;
            }
        }

        return q;
    }
};

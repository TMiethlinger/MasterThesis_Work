#include "ParticleAssignment.hpp"

#include <iostream>
#include <chrono>
#include <thread>
#include <sys/stat.h>

constexpr double ParticleAssignment::precision;
constexpr double ParticleAssignment::cost_max;

// Assignment methods
ParticleAssignment::Result ParticleAssignment::MinCostMatching(std::vector<std::vector<double>> &Ri, std::vector<std::vector<double>> &Rj, int M, int B)
{
    auto start = std::chrono::high_resolution_clock::now();
    int N = Ri.size();

    if(N != Rj.size())
    {
        // No assignment possible
        std::vector<int> q_final(N, -1);
        std::vector<int> q_inv_final(N, -1);
        std::vector<double> cost_vector;
        ParticleAssignment::CostObject cost(ParticleAssignment::cost_max, cost_vector);
        std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;
        ParticleAssignment::Result res(false, cost, q_final, q_inv_final, elapsed.count());
        return res;
    }
    std::vector<int> q_final(N, -1);
    std::vector<int> q_inv_final(N, -1);

    // Store all results in the list
    std::list<ParticleAssignment::Result> partialresults;

    // Create costskeleton
    std::vector<std::vector<ParticleAssignment::SkeletonNode>> costskeleton = ParticleAssignment::CreateCostSkeleton(Ri, Rj, M);
    // Assign CostSkeleton using a minimum cost metric, and add to resultslist
    ParticleAssignment::Result result = ParticleAssignment::MinCostMatchingSkeleton(costskeleton);
    partialresults.push_back(result);
	
    // Count number of unassigned nodes
    int countunmatched = std::count_if(result.q.begin(), result.q.end(), [](int i){ return i == -1; });

    std::vector<std::vector<double>> Ri_prev = Ri;
    std::vector<std::vector<double>> Rj_prev = Rj;

    while(countunmatched > B)
    {
        std::vector<int> q = result.q;
        std::vector<int> q_inv = result.q_inv;

        std::vector<std::vector<double>> Ri_mod(countunmatched);
        std::vector<std::vector<double>> Rj_mod(countunmatched);

        int N_prev = Ri_prev.size();
        for(int i = 0, j = 0, k = 0; i < N_prev; i++)
        {
            if(q[i] == -1)
            {
                Ri_mod[j] = Ri_prev[i];
                j++;
            }
            if(q_inv[i] == -1)
            {
                Rj_mod[k] = Rj_prev[i];
                k++;
            }
        }

        costskeleton = ParticleAssignment::CreateCostSkeleton(Ri_mod, Rj_mod, M);

        // Assign CostSkeleton using a minimum cost metric, and add to resultslist
        result = ParticleAssignment::MinCostMatchingSkeleton(costskeleton);
        partialresults.push_back(result);
        countunmatched = std::count_if(result.q.begin(), result.q.end(), [](int i){ return i == -1; });
        Ri_prev = Ri_mod;
        Rj_prev = Rj_mod;
    }

    if(countunmatched != 0) // we know: countunmatched <= B
    {
        std::vector<int> q = result.q;
        std::vector<int> q_inv = result.q_inv;

        std::vector<std::vector<double>> Ri_mod(countunmatched);
        std::vector<std::vector<double>> Rj_mod(countunmatched);

        int N_prev = Ri_prev.size();
        for(int i = 0, j = 0, k = 0; i < N_prev; i++)
        {
            if(q[i] == -1)
            {
                Ri_mod[j] = Ri_prev[i];
                j++;
            }
            if(q_inv[i] == -1)
            {
                Rj_mod[k] = Rj_prev[i];
                k++;
            }
        }

        std::vector<std::vector<double>> costmatrix = ParticleAssignment::CreateCostMatrix(Ri_mod, Rj_mod);

        // Assign CostSkeleton using a minimum cost metric, and add to resultslist
        result = ParticleAssignment::MinCostMatchingMatrix(costmatrix);
        partialresults.push_back(result);
    }

    for(std::list<ParticleAssignment::Result>::iterator it = partialresults.begin(); it != partialresults.end(); it++)
    {
        for(int i = 0, j = 0, k = 0; i < N && j < (*it).q.size(); i++)
        {
            if(q_final[i] == -1)
            {
                q_final[i] = (*it).q[j];
                j++;
            }
            if(q_inv_final[i] == -1)
            {
                q_inv_final[i] = (*it).q[k];
                k++;
            }
        }
    }

    ParticleAssignment::CostObject cost = ParticleAssignment::CostFromMatchingVectors(Ri, Rj, q_final);
    bool success = cost.total_cost < ParticleAssignment::cost_max ? true : false;
    std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;
    ParticleAssignment::Result res(success, cost, q_final, q_inv_final, elapsed.count());

    return res;
}
ParticleAssignment::Result ParticleAssignment::MinCostMatchingMatrix(std::vector<std::vector<double>> &costmatrix)
{
    auto start = std::chrono::high_resolution_clock::now();
    int N = costmatrix.size();
    if(N != costmatrix[0].size())
    {
        std::vector<int> q = std::vector<int>(N, -1);
        std::vector<int> q_inv = std::vector<int>(N, -1);
        std::vector<double> cost_vector;
        ParticleAssignment::CostObject cost(ParticleAssignment::cost_max, cost_vector);
        std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;
        ParticleAssignment::Result res(false, cost, q, q_inv, elapsed.count());
        return res;
    }

    // construct dual feasible solution
    std::vector<double> u(N, ParticleAssignment::cost_max);
    std::vector<double> v(N, ParticleAssignment::cost_max);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            u[i] = std::min(u[i], costmatrix[i][j]);
        }
    }
    for (int j = 0; j < N; j++)
    {
        for (int i = 0; i < N; i++)
        {
            v[j] = std::min(v[j], costmatrix[i][j] - u[i]);
        }
    }

    // construct primal solution satisfying complementary slackness
    std::vector<int> q = std::vector<int>(N, -1);
    std::vector<int> q_inv = std::vector<int>(N, -1);
    int mated = 0;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (q_inv[j] != -1)
                continue;
            if (fabs(costmatrix[i][j] - u[i] - v[j]) < ParticleAssignment::precision)
            {
                q[i] = j;
                q_inv[j] = i;
                mated++;
                break;
            }
        }
    }

    std::vector<double> dist(N);
    std::vector<int> dad(N);
    std::vector<int> seen(N);

    // repeat until primal solution is feasible
    while (mated < N)
    {
        // find an unmatched left node
        int s = 0;
        while (q[s] != -1)
        {
            s++;
        }

        // initialize Dijkstra
        std::fill(dad.begin(), dad.end(), -1);
        std::fill(seen.begin(), seen.end(), 0);
        for (int k = 0; k < N; k++)
        {
            dist[k] = costmatrix[s][k] - u[s] - v[k];
        }

        int j = 0;
        while (true)
        {
            // find closest
            j = -1;
            for (int k = 0; k < N; k++)
            {
                if (!seen[k])
                {
                    if (j == -1 || dist[k] < dist[j])
                    {
                        j = k;
                    }
                }
            }
            seen[j] = 1;

            // termination condition
            if (q_inv[j] == -1)
                break;

            // relax neighbours
            const int i = q_inv[j];
            for (int k = 0; k < N; k++)
            {
                if (!seen[k])
                {
                    const double new_dist = dist[j] + costmatrix[i][k] - u[i] - v[k];
                    if (new_dist < dist[k])
                    {
                        dist[k] = new_dist;
                        dad[k] = j;
                    }
                }
            }
        }

        // update dual variables
        for (int k = 0; k < N; k++)
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

    ParticleAssignment::CostObject cost = ParticleAssignment::CostFromMatchingMatrix(costmatrix, q);
    bool success = cost.total_cost < ParticleAssignment::cost_max ? true : false;
    std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;
    ParticleAssignment::Result res(success, cost, q, q_inv, elapsed.count());

    return res;
}
ParticleAssignment::Result ParticleAssignment::MinCostMatchingSkeleton(std::vector<std::vector<ParticleAssignment::SkeletonNode>> &costskeleton)
{
    auto start = std::chrono::high_resolution_clock::now();
    int N = costskeleton.size();

    std::function<bool(std::pair<int, double>, std::pair<int, double>)> compIntDoublePairSet =
        [](std::pair<int, double> x1, std::pair<int, double> x2)
        {
            if(x1.second < x2.second)
            {
                return true;
            }
            else if(x1.second > x2.second)
            {
                return false;
            }
            else
            {
                return x1.first < x2.first;
            }
        };

    // construct dual feasible solution
    std::vector<double> u(N, ParticleAssignment::cost_max);
    std::vector<double> v(N, ParticleAssignment::cost_max);
    for (int i = 0; i < N; i++)
    {
        const int size = costskeleton[i].size();
        for (int j = 0; j < size; j++)
        {
            u[i] = std::min(u[i], costskeleton[i][j].cost);
        }
    }
    for (int i = 0; i < N; i++)
    {
        const int size = costskeleton[i].size();
        for (int j = 0, y; j < size; j++)
        {
            y = costskeleton[i][j].index_column;
            v[y] = std::min(v[y], costskeleton[i][j].cost - u[i]);
        }
    }

    // construct primal solution satisfying complementary slackness
    std::vector<int> q = std::vector<int>(N, -1);
    std::vector<int> q_inv = std::vector<int>(N, -1);
    int mated = 0;
    for (int i = 0; i < N; i++)
    {
        const int size = costskeleton[i].size();
        for (int j = 0, y; j < size; j++)
        {
            y = costskeleton[i][j].index_column;
            if (q_inv[y] == -1)
            {
                if (fabs(costskeleton[i][j].cost - u[i] - v[y]) < ParticleAssignment::precision)
                {
                    q[i] = y;
                    q_inv[y] = i;
                    mated++;
                    break;
                }
            }
        }
    }

    std::vector<double> dist(N);
    std::vector<int> dad(N);
    std::vector<int> seen(N);

    // Memory variable for s: allows to skip non-matched nodes,
    // and go to the next one
    std::list<int> slist;
    for(int i = 0; i < N; i++)
    {
        slist.push_back(i);
    }

    // repeat until primal solution is feasible
    // for (int ctr1 = 0; mated < n; ctr1++)
    int s = -1;
    std::list<int>::iterator iter_s;
    std::set<std::pair<int, double>, std::function<bool(std::pair<int, double>, std::pair<int, double>)>> id_set(compIntDoublePairSet);
    std::pair<bool, std::set<std::pair<int, double>, std::function<bool(std::pair<int, double>, std::pair<int, double>)>>::iterator> default_pair = std::make_pair(false, id_set.end());
    std::vector<std::pair<bool, std::set<std::pair<int, double>, std::function<bool(std::pair<int, double>, std::pair<int, double>)>>::iterator>> bi_vector(N, default_pair);
    // repeat until primal solution is feasible
    while (mated < N)
    {
        // if: there was an assignment previously and consequently s == -1
        if(s == -1)
        {
            // find an unmatched left node
            iter_s = slist.begin();
            while(q[*iter_s] != -1 && iter_s != slist.end())
            {
                iter_s++;
            }
            if(iter_s == slist.end())
            {
                // no other assignments free
                std::vector<double> cost_vector;
                ParticleAssignment::CostObject cost(ParticleAssignment::cost_max, cost_vector);
                std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;
                ParticleAssignment::Result res(false, cost, q, q_inv, elapsed.count());
                return res;
            }
            else
            {
                s = *iter_s;
            }
        }
        // else: Use previously found unmatched left node

        // initialize Dijkstra
        std::fill(dad.begin(), dad.end(), -1);
        std::fill(seen.begin(), seen.end(), 0);
        std::fill(dist.begin(), dist.end(), ParticleAssignment::cost_max);
        const int size = costskeleton[s].size();

        for (int l = 0, y; l < size; l++)
        {
            y = costskeleton[s][l].index_column;
            dist[y] = costskeleton[s][l].cost - u[s] - v[y];
            if(bi_vector[y].first)
            {
                // The index y is already active:
                // 1. Delete existing entry
                id_set.erase(bi_vector[y].second);
                // 2. Add entry and set iterator
                bi_vector[y].second = id_set.insert(std::make_pair(y, dist[y])).first;
            }
            else
            {
                // The index y is inactive
                // 1. Set index active
                bi_vector[y].first = true;
                // 2. Add entry and set iterator
                bi_vector[y].second = id_set.insert(std::make_pair(y, dist[y])).first;
            }
        }

        int j;
        while (true)
        {
            j = -1;
            std::set<std::pair<int, double>>::iterator iter_min = id_set.begin();
            if(iter_min != id_set.end())
            {
                j = iter_min->first;
                bi_vector[j].first = false;
                bi_vector[j].second = id_set.end();
                id_set.erase(iter_min);
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
                    std::vector<double> cost_vector;
                    ParticleAssignment::CostObject cost(ParticleAssignment::cost_max, cost_vector);
                    std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;
                    ParticleAssignment::Result res(false, cost, q, q_inv, elapsed.count());
                    return res;
                }
                else
                {
                    s = *iter_s;
                    break;
                }
            }

            // termination condition
            if (q_inv[j] == -1)
                break;

            const int i = q_inv[j];
            const int size = costskeleton[i].size();
            for (int l = 0, y; l < size; l++)
            {
                y = costskeleton[i][l].index_column;
                if (!seen[y])
                {
                    const double new_dist = dist[j] + costskeleton[i][l].cost - u[i] - v[y];
                    if (new_dist < dist[y])
                    {
                        dist[y] = new_dist;
                        dad[y] = j;
                        if(bi_vector[y].first)
                        {
                            // The index y is already active:
                            // 1. Delete existing entry
                            id_set.erase(bi_vector[y].second);
                            // 2. Add entry and set iterator
                            bi_vector[y].second = id_set.insert(std::make_pair(y, dist[y])).first;
                        }
                        else
                        {
                            // The index y is inactive
                            // 1. Set index active
                            bi_vector[y].first = true;
                            // 2. Add entry and set iterator
                            bi_vector[y].second = id_set.insert(std::make_pair(y, dist[y])).first;
                        }
                    }
                }
            }
        }

        if(j != -1)
        {
            // update dual variables
            for (int k = 0; k < N; k++)
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
            id_set.clear();
            std::fill(bi_vector.begin(), bi_vector.end(), default_pair);
            s = -1;
        }
    }
    ParticleAssignment::CostObject cost = ParticleAssignment::CostFromMatchingSkeleton(costskeleton, q);
    std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;
    ParticleAssignment::Result res(true, cost, q, q_inv, elapsed.count());
    return res;
}

// Creator methods
std::vector<std::vector<double>> ParticleAssignment::CreateCostMatrix(std::vector<std::vector<double>> &Ri, std::vector<std::vector<double>> &Rj)
{
    int N = Ri.size();
    std::vector<std::vector<double>> costmatrix(N, std::vector<double>(N));

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            costmatrix[i][j] = EuclideanMetric(Ri[i], Rj[j]);
        }
    }

    return costmatrix;
}
std::vector<std::vector<double>> ParticleAssignment::CreateCostMatrix(std::vector<std::vector<double>> &Ri, std::vector<std::vector<double>> &Rj, int M)
{
    int N = Ri.size();
    std::vector<std::vector<double>> costmatrix(N, std::vector<double>(N));

    std::vector<std::pair<double, int>> pairskeleton_row(N);
    std::vector<std::pair<double, int>> pairskeleton_column(N);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            std::pair<double, int> p;
            p.first = EuclideanMetric(Ri[i], Rj[j]);
            p.second = j;
            pairskeleton_row[j] = p;
            costmatrix[i][j] = ParticleAssignment::cost_max;
        }

        std::nth_element(pairskeleton_row.begin(), pairskeleton_row.begin() + M, pairskeleton_row.end(), [](std::pair<double, int> a, std::pair<double, int> b) {return a.first < b.first;});

        for (int j = 0; j < M; j++)
        {
            costmatrix[i][pairskeleton_row[j].second] = pairskeleton_row[j].first;
        }
    }

    for (int j = 0; j < N; j++)
    {
        for (int i = 0; i < N; i++)
        {
            std::pair<double, int> p;
            p.first = EuclideanMetric(Ri[i], Rj[j]);
            p.second = i;
            pairskeleton_column[i] = p;
        }

        std::nth_element(pairskeleton_column.begin(), pairskeleton_column.begin() + M, pairskeleton_column.end(), [](std::pair<double, int> a, std::pair<double, int> b) {return a.first < b.first;});

        for (int i = 0; i < M; i++)
        {
            costmatrix[pairskeleton_column[i].second][j] = pairskeleton_column[i].first;
        }
    }

    return costmatrix;
}
std::vector<std::vector<ParticleAssignment::SkeletonNode>> ParticleAssignment::CreateCostSkeleton(std::vector<std::vector<double>> &Ri, std::vector<std::vector<double>> &Rj, int M)
{
    int N = Ri.size();
    std::vector<std::vector<ParticleAssignment::SkeletonNode>> costskeleton(N);
    // Find the M nearest neighbours per row
    std::vector<std::vector<std::pair<double, int>>> costskeleton_row_only(N, std::vector<std::pair<double, int>>(M));
    std::vector<std::pair<double, int>> pairskeleton_row(N);
    std::vector<std::pair<double, int>> pairskeleton_row_part(M);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            pairskeleton_row[j] = std::make_pair(EuclideanMetric(Ri[i], Rj[j]), j);
        }

        std::nth_element(pairskeleton_row.begin(), pairskeleton_row.begin() + M, pairskeleton_row.end(), [](std::pair<double, int> a, std::pair<double, int> b) {return a.first < b.first;});
        std::copy(pairskeleton_row.begin(), pairskeleton_row.begin() + M, pairskeleton_row_part.begin());
        std::sort(pairskeleton_row_part.begin(), pairskeleton_row_part.end(), [](std::pair<double, int> a, std::pair<double, int> b) {return a.second < b.second;});
        std::copy(pairskeleton_row_part.begin(), pairskeleton_row_part.end(), costskeleton_row_only[i].begin());
    }

    // Find the M nearest neighbours per column
    std::vector<std::vector<std::pair<double, int>>> costskeleton_column_only(N, std::vector<std::pair<double, int>>(M));
    std::vector<std::pair<double, int>> pairskeleton_column(N);
    std::vector<std::pair<double, int>> pairskeleton_column_part(M);
    for (int j = 0; j < N; j++)
    {
        for (int i = 0; i < N; i++)
        {
            pairskeleton_column[i] = std::make_pair(EuclideanMetric(Ri[i], Rj[j]), i);
        }

        std::nth_element(pairskeleton_column.begin(), pairskeleton_column.begin() + M, pairskeleton_column.end(), [](std::pair<double, int> a, std::pair<double, int> b) {return a.first < b.first;});
        std::copy(pairskeleton_column.begin(), pairskeleton_column.begin() + M, pairskeleton_column_part.begin());
        std::sort(pairskeleton_column_part.begin(), pairskeleton_column_part.end(), [](std::pair<double, int> a, std::pair<double, int> b) {return a.second < b.second;});
        std::copy(pairskeleton_column_part.begin(), pairskeleton_column_part.end(), costskeleton_column_only[j].begin());
    }

    // Merge such that the number of neighbours is at least m row-&columnwise
    std::vector<std::list<std::pair<double, int>>> costskeleton_list_full(N, std::list<std::pair<double, int>>());
    for (int i = 0; i < N; i++)
    {
        costskeleton_list_full[i].assign(costskeleton_row_only[i].begin(), costskeleton_row_only[i].end());
    }
    int index_row;
    bool found;
    for (int j = 0; j < N; j++)
    {
        for(std::vector<std::pair<double, int>>::iterator it = costskeleton_column_only[j].begin(); it != costskeleton_column_only[j].end(); it++)
        {
            index_row = it->second;
            found = false;
            for(std::vector<std::pair<double, int>>::iterator jt = costskeleton_row_only[index_row].begin(); jt != costskeleton_row_only[index_row].end(); jt++)
            {
                if(j == jt->second)
                {
                    found = true;
                    break;
                }
            }
            if(!found)
            {
                costskeleton_list_full[index_row].push_back(std::make_pair(it->first, j));
            }
        }
    }

    // Sort the merged lists
    for(int i = 0; i < N; i++)
    {
        costskeleton_list_full[i].sort([](std::pair<double, int> a, std::pair<double, int> b) {return a.second < b.second;});
    }
    int size;

    // Set the skeletons to the computed neighbours
    for(int i = 0; i < N; i++)
    {
        size = costskeleton_list_full[i].size();
        costskeleton[i].resize(size);
        {
            int j = 0;
            for(std::list<std::pair<double, int>>::iterator jt = costskeleton_list_full[i].begin(); jt != costskeleton_list_full[i].end(); jt++, j++)
            {
                ParticleAssignment::SkeletonNode s(jt->first, jt->second);
                costskeleton[i][j] = s;
            }
        }
    }

    return costskeleton;
}
/* std::vector<std::vector<double>> ParticleAssignment::CreateCostMatrixPBC(std::vector<std::vector<double>> &Ri, std::vector<std::vector<double>> &Rj, double L, double L2)
{
    int N = Ri.size();
    std::vector<std::vector<double>> costmatrix(N, std::vector<double>(N));

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            costmatrix[i][j] = EuclideanMetricPBC(Ri[i], Rj[j], L, L2);
        }
    }

    return costmatrix;
} */

// Compute cost
ParticleAssignment::CostObject ParticleAssignment::CostFromMatchingVectors(std::vector<std::vector<double>> &Ri, std::vector<std::vector<double>> &Rj, std::vector<int> &q)
{
    int N = Ri.size();

    double total_cost = 0;
    std::vector<double> cost_vector(N);

    double d;
    for (int i = 0; i < N; i++)
    {
        d = EuclideanMetric(Ri[i], Rj[q[i]]);
        total_cost += d;
        cost_vector[i] = d;
    }

    ParticleAssignment::CostObject cost(total_cost, cost_vector);
    return cost;
}
ParticleAssignment::CostObject ParticleAssignment::CostFromMatchingMatrix(std::vector<std::vector<double>> &costmatrix, std::vector<int> &q)
{
    int N = costmatrix.size();

    double total_cost = 0;
    std::vector<double> cost_vector(N);

    double d;
    for (int i = 0; i < N; i++)
    {
        d = costmatrix[i][q[i]];
        total_cost += costmatrix[i][q[i]];
        cost_vector[i] = d;
    }

    ParticleAssignment::CostObject cost(total_cost, cost_vector);
    return cost;
}
ParticleAssignment::CostObject ParticleAssignment::CostFromMatchingSkeleton(std::vector<std::vector<ParticleAssignment::SkeletonNode>> &costskeleton, std::vector<int> &q)
{
    int N = costskeleton.size(), M;

    double total_cost = 0;
    std::vector<double> cost_vector(N);
    bool found;

    for (int i = 0; i < N; i++)
    {
        found = false;
        M = costskeleton[i].size();
        for (int j = 0; j < M; j++)
        {
            if (q[i] == costskeleton[i][j].index_column)
            {
                total_cost += costskeleton[i][j].cost;
                cost_vector[i] = costskeleton[i][j].cost;
                found = true;
                break;
            }
        }
        if(!found)
        {
            total_cost += ParticleAssignment::cost_max;
            break;
        }
    }

    ParticleAssignment::CostObject cost(total_cost, cost_vector);
    return cost;
}

// Metric

double ParticleAssignment::EuclideanMetric(std::vector<double> &ri, std::vector<double> &rj)
{
    double dx, dy, dz;
    dx = fabs(ri[0] - rj[0]);
    dy = fabs(ri[1] - rj[1]);
    dz = fabs(ri[2] - rj[2]);
    return sqrt(dx*dx + dy*dy + dz*dz);
}


/*double ParticleAssignment::EuclideanMetric(std::vector<double> &ri, std::vector<double> &rj)
{
    double d;
	double sum = 0.0;
	for(int i = 0; i < ri.size(); i++)
	{
		d = ri[i] - rj[i];
		sum += d*d;
	}
    return sqrt(sum);
}
double ParticleAssignment::EuclideanMetricPBC(std::vector<double> &ri, std::vector<double> &rj, double L, double L2)
{
    double d;
	double sum = 0.0;
	for(int i = 0; i < ri.size(); i++)
	{
		d = fabs(ri[i] - rj[i]);
		d = d > L2 ? L - d : d;
		sum += d*d;
	}
    return sqrt(sum);
}*/

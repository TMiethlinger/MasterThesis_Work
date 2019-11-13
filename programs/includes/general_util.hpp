// Author:     Thomas Miethlinger B.Sc.
// Date:       16.10.2019
// Email:      t.miethlinger@protonmail.com
// University: Johannes Kepler University Linz, Austria

// STL
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <thread>
#include <sys/stat.h>

// Boost
#include <boost/algorithm/string.hpp>

#pragma once

namespace general_util
{
    std::vector<int> create_jobs_vector(int njobs, int world_rank, int world_size)
    {
        std::vector<int> job_vector;
        int njobs_per_rank_min = njobs / world_size;
        int njobs_per_rank_max = njobs % world_size == 0 ? njobs_per_rank_min : njobs_per_rank_min + 1;
        int min_rank = njobs % world_size;
        int njobs_my_rank = world_rank < min_rank ? njobs_per_rank_max : njobs_per_rank_min;
        job_vector.resize(njobs_my_rank);

        int ijobmin, ijobmax;
        if(world_rank < min_rank)
        {
            ijobmin = world_rank * njobs_per_rank_max;
            ijobmax = ijobmin + njobs_per_rank_max - 1;
        }
        else
        {
            ijobmin = min_rank * njobs_per_rank_max + (world_rank - min_rank) * njobs_per_rank_min;
            ijobmax = ijobmin + njobs_per_rank_min - 1;
        }

        for(int i = ijobmin, j = 0; i <= ijobmax; i++, j++)
        {
            job_vector[j] = i;
        }

        return job_vector;
    }

    std::vector<std::vector<int>> partition_total_jobs_vector(int njobs, int world_size)
    {
        std::vector<std::vector<int>> job_vector_partitioned(world_size);
        int njobs_per_rank_min = njobs / world_size;
        int njobs_per_rank_max = njobs % world_size == 0 ? njobs_per_rank_min : njobs_per_rank_min + 1;
        int min_rank = njobs % world_size;

        for(int r = 0; r < world_size; r++)
        {
            int njobs_my_rank = r < min_rank ? njobs_per_rank_max : njobs_per_rank_min;
            job_vector_partitioned[r].resize(njobs_my_rank);
            int ijobmin, ijobmax;

            if(r < min_rank)
            {
                ijobmin = r * njobs_per_rank_max;
                ijobmax = ijobmin + njobs_per_rank_max - 1;
            }
            else
            {
                ijobmin = min_rank * njobs_per_rank_max + (r - min_rank) * njobs_per_rank_min;
                ijobmax = ijobmin + njobs_per_rank_min - 1;
            }

            for(int i = ijobmin, j = 0; i <= ijobmax; i++, j++)
            {
                job_vector_partitioned[r][j] = i;
            }
        }

        return job_vector_partitioned;
    }

    std::vector<std::vector<int>> partition_total_jobs_vector(std::vector<int> job_vector, int world_size)
    {
        int njobs = job_vector.size();
        std::vector<std::vector<int>> job_vector_partitioned(world_size);
        int njobs_per_rank_min = njobs / world_size;
        int njobs_per_rank_max = njobs % world_size == 0 ? njobs_per_rank_min : njobs_per_rank_min + 1;
        int min_rank = njobs % world_size;

        for(int r = 0; r < world_size; r++)
        {
            int njobs_my_rank = r < min_rank ? njobs_per_rank_max : njobs_per_rank_min;
            job_vector_partitioned[r].resize(njobs_my_rank);
            int ijobmin, ijobmax;

            if(r < min_rank)
            {
                ijobmin = r * njobs_per_rank_max;
                ijobmax = ijobmin + njobs_per_rank_max - 1;
            }
            else
            {
                ijobmin = min_rank * njobs_per_rank_max + (r - min_rank) * njobs_per_rank_min;
                ijobmax = ijobmin + njobs_per_rank_min - 1;
            }

            for(int i = ijobmin, j = 0; i <= ijobmax; i++, j++)
            {
                job_vector_partitioned[r][j] = job_vector[i];
            }
        }

        return job_vector_partitioned;
    }


    int job_index_from_time_indices(int i, int j, int n)
    {
        int index = 0;
        for(int k = 0; k < i; k++)
        {
            index += (n - k - 1);
        }
        index += j - (i + 1);
        return index;
    }
};

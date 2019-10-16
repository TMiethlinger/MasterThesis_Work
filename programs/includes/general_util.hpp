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
    // MPI functions
 
    // Evenly assign jobs per rank
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
};

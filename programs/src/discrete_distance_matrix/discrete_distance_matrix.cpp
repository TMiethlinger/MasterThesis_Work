// Author:     Thomas Miethlinger B.Sc.
// Date:       29.10.2019
// Email:      t.miethlinger@protonmail.com
// University: Johannes Kepler University Linz, Austria

// STL
#include <algorithm> // std::fill
#include <array>
#include <cmath> // sqrt
#include <iostream>
#include <string>
#include <utility> // std::pair
#include <vector>

// MPI
#include <mpi.h>

// Boost
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

// Internal classes
#include "../../includes/io_util.hpp"
#include "../../includes/general_util.hpp"
#include "../hungarian_algorithm/ha_core.hpp"
#include "../hungarian_algorithm/ha_costobject.hpp"
#include "../hungarian_algorithm/ha_distance.hpp"

using std::array;
using std::pair;
using std::size_t;
using std::string;
using std::vector;
using namespace boost::program_options;

typedef vector<int> VI;
typedef vector<double> VD;
typedef vector<VD> VVD;
typedef array<int, 2> AI2;
typedef pair<int, double> PID;

// Variables for the processing of the input dataset
// They are set by the programm arguments argv
int dim; // Dimension of a single-particle mechanical state
int N; // Number of particles
int N_match; // Number of particles to consider for matching
double l_inv; // Inverse of characteristic lengthscale
double v_inv; // Inverse of characteristic velocityscale
int tmin;
int tmax;
int tstep;
int nsteps;
int njobs;

// IO
string inputfolder_parent = "/home/k3501/k354524/master_thesis_work/data/";
string inputfolder_relative; // e.g. Liggghts/.../
string inputfilenamepart;
string outputfolder_parent = "/home/k3501/k354524/master_thesis_work/results/field_distance_matrix/";
string outputfolder_relative; // e.g. Liggghts/.../
string outputfilename = "field_distance_matrix.txt";

int main(int argc, char * argv[])
{
    // Initialize MPI environment
    MPI_Init(&argc, &argv);

    // Query MPI environment for rank and size
    int world_size;
    int world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Parse command line using boost::program_options
    options_description desc_commandline;
    desc_commandline.add_options()
    ("dim", value<int>()->default_value(6), "Dimension of a single particle state")    
    ("N", value<int>()->default_value(1), "Particle number N")
    ("N_match", value<int>()->default_value(1), "Number of particle for matching N_match")
    ("l_inv", value<double>()->default_value(1.0), "Inverse of characteristic lengthscale")
    ("v_inv", value<double>()->default_value(1.0), "Inverse of characteristic velocityscale")
    ("tmin", value<int>()->default_value(0), "First time index")
    ("tmax", value<int>()->default_value(1), "Last time index")
    ("tstep", value<int>()->default_value(1), "Time index step size")
    ("inputfolder_relative", value<string>()->default_value("Liggghts/"), "Relative input foldername")
    ("inputfilenamepart", value<string>()->default_value("dump"), "Input file name part")
    ("outputfolder_relative", value<string>()->default_value("Liggghts/"), "Relative output foldername");
    variables_map vm;
    store(parse_command_line(argc, argv, desc_commandline), vm);
    N = vm["N"].as<int>();
    N_match = vm["N_match"].as<int>();
    l_inv = vm["l_inv"].as<double>();
    v_inv = vm["v_inv"].as<double>();
    tmin = vm["tmin"].as<int>();
    tmax = vm["tmax"].as<int>();
    tstep = vm["tstep"].as<int>();
    nsteps = (tmax - tmin) / tstep + 1;
    njobs = (nsteps*nsteps - nsteps)/2;
    inputfolder_relative = vm["inputfolder_relative"].as<string>();
    inputfilenamepart = vm["inputfilenamepart"].as<string>();
    outputfolder_relative = vm["outputfolder_relative"].as<string>();

    // Let rank 0 create the outputdirectory, if it does not exist
    if(world_rank == 0)
    {
        string outputfolder = outputfolder_parent + outputfolder_relative;
        io_util::create_working_directory(outputfolder.c_str());
    }

    // Get all jobs, i.e. int i which produce a pair of {i_1, i_2} for which the distance d_{i_1,i_2} is to be calculated
    VI jobs_vector_my_rank = general_util::create_jobs_vector(njobs, world_rank, world_size);
    size_t njobs_my_rank = jobs_vector_my_rank.size();
    // Define a vector total_jobs_vector which transforms from i to {t_1, t_2} for all jobs i
    vector<AI2> total_jobs_vector(njobs);

    for(int t1 = tmin, i = 0; t1 <= (tmax - tstep); t1 += tstep)
    {
        for(int t2 = t1 + tstep; t2 <= tmax; t2 += tstep, i++)
        {
            total_jobs_vector[i] = {t1, t2};
        }
    }

    // For each rank and its jobs, save the result of the computation into field_distance_results
    VD discrete_distance_results(njobs_my_rank, 0);

    // Initialize particle state vectors Xa and Xb
    VVD Xa(N, VD(dim));
    VVD Xb(N, VD(dim));

    // Variables for the indices of each job
    int t1, t2, t1_prev, t2_prev;
    t1 = t2 = t1_prev = t2_prev = -1;
    array<int, 2> job_index_arr;

    for(size_t i = 0; i != njobs_my_rank; i++)
    {
        job_index_arr = total_jobs_vector[jobs_vector_my_rank[i]];

        // Files may be saved per time-step with an offset tmin and scaling tstep.
        t1 = job_index_arr[0];
        // In order to avoid two-times reading from the same file
        if(t1 != t1_prev)
        {
            // Read input dump file
            if (!io_util::read_positions_3D(inputfolder_parent + inputfolder_relative + inputfilenamepart + std::to_string(t1), N, Xa))
            {
                std::cout << "Error: Could not read file " << (inputfolder_parent + inputfolder_relative + inputfilenamepart + std::to_string(t1)) << "!" << std::endl;
                return 1;
            }
        }

        // Same for state b
        t2 = job_index_arr[1];
        if(t2 != t2_prev)
        {
            if (!io_util::read_positions_3D(inputfolder_parent + inputfolder_relative + inputfilenamepart + std::to_string(t2), N, Xb))
            {
                std::cout << "Error: Could not read file " << (inputfolder_parent + inputfolder_relative + inputfilenamepart + std::to_string(t1)) << "!" << std::endl;
                return 1;
            }
        }

        // Compute the distance adjacency lists
        vector<vector<PID>> cost_adjlist = ha_costobject::create_costobject_adjlist_plain(N, N_match, Xa, Xb, ha_distance::d_sum_3, l, v);

        // Compute the hungarian algorithm matching solution for adjacency lists
        VI q = ha_core::min_cost_matching_adjlist(cost_adjlist);

        // Compute the total cost from the matching solution
        discrete_distance_results[i] = compute_cost_from_adjlist(N, cost_adjlist, q);

        // Update indices
        t1_prev = t1;
        t2_prev = t2;
    }

    // Write the intermediate distance results 
    // IO::WriteVectorResult(outputfolder_parent + outputfolder_relative + "field_distance_matrix_" + std::to_string(world_rank) + ".txt", field_distance_results);

    // Rank 0 needs to compose the field_distance_matrix
    if(world_rank == 0)
    {
        // Save all distances from each rank into the following vector
        VD total_field_distance_results(njobs, 0);
        // Calculate how the full problem (computing the distance_matrix) was split up into jobs
        int njobs_per_rank_min = njobs / world_size;
        int njobs_per_rank_max = njobs % world_size == 0 ? njobs_per_rank_min : njobs_per_rank_min + 1;
        int min_rank = njobs % world_size;

        // Add own intermediate results to total resul vector
        for(int i = 0; i < njobs_per_rank_max; i++)
        {
            total_field_distance_results[i] = field_distance_results[i];
        }

        // Receive all intermediate results, with vector lengths recv_result_size.
        // recv_offset is a variable which shifts the beginning of the receive buffer.
        int recv_result_size;
        int recv_offset = njobs_per_rank_max;
        for(int w = 1; w < world_size; w++)
        {
            recv_result_size = w < min_rank ? njobs_per_rank_max : njobs_per_rank_min;
            MPI_Recv(total_field_distance_results.data() + recv_offset, recv_result_size, MPI_DOUBLE, w, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            recv_offset += recv_result_size;
        }

        // Set the final field_distance_matrix from received results
        VVD field_distance_matrix(nsteps, VD(nsteps, 0));
        for(int i1 = 0, i = 0; i1 < nsteps - 1; i1++)
        {
            for(int i2 = i1 + 1; i2 < nsteps; i2++, i++)
            {
                field_distance_matrix[i1][i2] = field_distance_matrix[i2][i1] = total_field_distance_results[i];
            }
        }

        // Write the computed distance matrix based on fields into this filepathname
        io_util::write_matrix_result(outputfolder_parent + outputfolder_relative + outputfilename, field_distance_matrix);
    }
    // Each other rank sends their intermediate results to rank 0
    else
    {
        MPI_Send(field_distance_results.data(), field_distance_results.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    // Finalize the MPI environment.
    MPI_Finalize();

    return 0;
}

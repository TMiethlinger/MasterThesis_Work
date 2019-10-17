// Author:     Thomas Miethlinger B.Sc.
// Date:       16.10.2019
// Email:      t.miethlinger@protonmail.com
// University: Johannes Kepler University Linz, Austria

// STL
#include <iostream>
#include <string>
#include <cmath> // sqrt

#include <array>
#include <vector>
#include <cctype>
#include <cstdlib>
#include <unistd.h>

// MPI
#include <mpi.h>

// Boost
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

// Internal classes
#include "../rqa_core/rqa_core.hpp"
#include "../../../includes/general_util.hpp"
#include "../../../includes/io_util.hpp"

using std::size_t;
using std::string;
using std::array;
using std::vector;
using namespace boost::program_options;
typedef vector<int> VI;
typedef vector<VI> VVI;
typedef vector<double> VD;
typedef vector<VD> VVD;


// Variables for the processing of the input dataset
// They are set by the programm arguments argv
int n;
int n2;
double epsmin;
double epsmax;
double deps;
int neps;
int njobs;

// IO
string inputfolder_parent = "/home/k3501/k354524/master_thesis_work/results/";
string inputfolder_relative; // e.g. Liggghts/.../
string inputfilename;
string outputfolder_parent = "/home/k3501/k354524/master_thesis_work/results/recurrence_quantitative_analysis/recurrence_rate/";
string outputfolder_relative; // e.g. Liggghts/.../
string outputfilename;

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
    ("n", value<int>()->default_value(1), "Number of time steps n")
    ("epsmin", value<double>()->default_value(0.0), "Minimum epsilon for recurrence threshold")
    ("epsmax", value<double>()->default_value(1.0), "Maximum epsilon for recurrence threshold")
    ("neps", value<int>()->default_value(8), "Number of different epsilon for evaluation")
    ("inputfolder_relative", value<string>()->default_value("Liggghts/"), "Relative input foldername")
    ("inputfilename", value<string>()->default_value("distance_matrix.txt"), "Input filename")
    ("outputfolder_relative", value<string>()->default_value("Liggghts/"), "Relative output foldername")
    ("outputfilename", value<string>()->default_value("recurrence_rates.txt/"), "Output filename");

    variables_map vm;
    store(parse_command_line(argc, argv, desc_commandline), vm);
    n = vm["n"].as<int>();
    n2 = n * n;
    epsmin = vm["epsmin"].as<double>();
    epsmax = vm["epsmax"].as<double>();
    neps = vm["neps"].as<int>();
    njobs = neps;
    deps = (epsmax - epsmin) / neps;
    inputfolder_relative = vm["inputfolder_relative"].as<string>();
    inputfilename = vm["inputfilename"].as<string>();
    outputfolder_relative = vm["outputfolder_relative"].as<string>();
    outputfilename = vm["outputfilename"].as<string>();

    // Let rank 0 create the outputdirectory, if it does not exist
    if(world_rank == 0)
    {
        string outputfolder = outputfolder_parent + outputfolder_relative;
        io_util::create_working_directory(outputfolder.c_str());
    }

    // Declare variable for normalized distance matrix and recurrence matrix on each rank
    // All matrix variables are saved in row-major mode
    VD norm_distance_matrix(n2, 0.0);
    VD recurrence_matrix(n2, 0.0);

    // Let rank 0
    if(world_rank == 0)
    {
        // Declare variable for distance matrix
        VD distance_matrix(n2, 0.0);

        // Read the distance_matrix from the input file into the variable
        io_util::read_double_matrix_row_major(inputfolder_parent + inputfolder_relative + inputfilename, distance_matrix, " ", 0, n, 0, n);

        // Compute the corresponding normalized distance matrix
        rqa_core::normalize_distance_matrix(distance_matrix, norm_distance_matrix);

        // Send the normalized distance matrix to each other rank
        for(int w = 1; w < world_size; w++)
        {
            MPI_Send(norm_distance_matrix.data(), n2, MPI_DOUBLE, w, 0, MPI_COMM_WORLD);
        }
    }
    else
    {
        // Receive normalized distance matrix from rank 0
        MPI_Recv(norm_distance_matrix.data(), n2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // Get all jobs i, i.e. compute the recurrence rate w.r.t. eps(i) = i * deps
    VI jobs_vector = general_util::create_jobs_vector(njobs, world_rank, world_size);
    size_t njobs_my_rank = jobs_vector.size();

    // For each rank and its jobs, save the result of the computation into field_distance_results
    VD recurrence_rate_results(njobs_my_rank, 0);

    // Compute the all recurrence rates based on jobs_vector
    double eps;
    double rr;
    for(size_t i = 0; i != njobs_my_rank; i++)
    {
        eps = jobs_vector[i] * deps + deps / 2.0;
        rqa_core::calculate_recurrence_matrix(norm_distance_matrix, eps, recurrence_matrix);
        rr = rqa_core::recurrence_rate(recurrence_matrix);
        recurrence_rate_results[i] = rr;
    }

    // Write the intermediate distance results 
    // io_util::write_vector_result(outputfolder_parent + outputfolder_relative + "recurrence_rate_" + std::to_string(world_rank) + ".txt", recurrence_rate_results);

    // Rank 0 needs to compose the recurrence_rate_results
    if(world_rank == 0)
    {
        // Save all recurrence rates from each rank into the following vector
        VD total_recurrence_rate_results(njobs, 0);
        //Calculate how the full problem (computing the recurrence rates w.r.t. eps) was split up into jobs
        int njobs_per_rank_min = njobs / world_size;
        int njobs_per_rank_max = njobs % world_size == 0 ? njobs_per_rank_min : njobs_per_rank_min + 1;
        int min_rank = njobs % world_size;

        // Add own intermediate results to total result vector
        for(int i = 0; i < njobs_per_rank_max; i++)
        {
            total_recurrence_rate_results[i] = recurrence_rate_results[i];
        }

        // Receive all intermediate results, with vector lengths recv_result_size.
        // recv_offset is a variable which shifts the beginning of the receive buffer.
        int recv_result_size;
        int recv_offset = njobs_per_rank_max;
        for(int w = 1; w < world_size; w++)
        {
            recv_result_size = w < min_rank ? njobs_per_rank_max : njobs_per_rank_min;
            MPI_Recv(total_recurrence_rate_results.data() + recv_offset, recv_result_size, MPI_DOUBLE, w, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            recv_offset += recv_result_size;
        }


        // Create a two-dimensional vector for file writing
        VVD output_result_matrix(njobs, VD(2, 0));
        for(int i = 0; i < neps; i++)
        {
            output_result_matrix[i][0] = deps * i + deps / 2.0;
            output_result_matrix[i][1] = total_recurrence_rate_results[i];
        }
        

        // Write the computed recurrence rates into this filepathname
        io_util::write_matrix_result(outputfolder_parent + outputfolder_relative + outputfilename, output_result_matrix);
    }
    // Each other rank sends their intermediate results to rank 0
    else
    {
        MPI_Send(recurrence_rate_results.data(), recurrence_rate_results.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    // Finalize the MPI environment.
    MPI_Finalize();

    return 0;
}



// STL
#include <iostream>
#include <string>
#include <cmath> // sqrt
#include <algorithm> // std::fill
#include <numeric> // std::accumulate
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
#include "../recurrence_quantitative_analysis/recurrence_quantitative_analysis.hpp"
#include "../../general.hpp"

using std::size_t;
using std::string;
using std::array;
using std::vector;
using namespace boost::program_options;
typedef vector<int> VI;
typedef vector<VI> VVI;
typedef vector<double> VD;
typedef vector<VD> VVD;

// General constants
const double precision = 0.000000000001;
const int digits = 3;

// Variables for the processing of the input dataset
// They are set by the programm arguments argv
int N;
double epsmin;
double epsmax;
double deps;
int neps;
int njobs;

// Calculation
double recurrence_rate(VD &recurrence_matrix);

// IO
string inputfolder_parent = "/home/k3501/k354524/master_thesis_work/results/";
string inputfolder_relative; // e.g. Liggghts/.../
string inputfilename;
string outputfolder_parent = "/home/k3501/k354524/master_thesis_work/results/recurrence_quantitative_analysis/recurrence_rate/";
string outputfolder_relative; // e.g. Liggghts/.../
string outputfilename = "field_distance_matrix.txt";

// Jobs
vector<VI> CreateJobsVector(int world_rank, int world_size);

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
    ("N", value<int>()->default_value(1), "Number of time steps N")
    ("epsmin", value<double>()->default_value(0.0), "Minimum epsilon for recurrence threshold")
    ("epsmax", value<double>()->default_value(1.0), "Maximum epsilon for recurrence threshold")
    ("neps", value<int>()->default_value(8), "Number of different epsilon for evaluation")
    ("inputfolder_relative", value<string>()->default_value("Liggghts/"), "Relative input foldername")
    ("inputfilenamepart", value<string>()->default_value("distance_matrix.txt"), "Input file name part")
    ("outputfolder_relative", value<string>()->default_value("Liggghts/"), "Relative output foldername");
    variables_map vm;
    store(parse_command_line(argc, argv, desc_commandline), vm);
    N = vm["N"].as<int>();
    epsmin = vm["epsmin"].as<double>();
    epsmax = vm["epsmax"].as<double>();
    neps = vm["neps"].as<int>();
    deps = (epsmax - epsmin) / neps;
    inputfolder_relative = vm["inputfolder_relative"].as<string>();
    inputfilenamepart = vm["inputfilenamepart"].as<string>();
    outputfolder_relative = vm["outputfolder_relative"].as<string>();

    // Let rank 0 create the outputdirectory, if it does not exist
    if(world_rank == 0)
    {
        string outputfolder = outputfolder_parent + outputfolder_relative;
        IO::CreateWorkingDirectory(outputfolder.c_str());
    }

    // Declare variable for normalized distance matrix and recurrence matrix on each rank
    VVD norm_distance_matrix(N, VD(N, 0));
    VVD recurrence_matrix(N, VD(N, 0));

    // Let rank 0
    if(world_rank == 0)
    {
        // Declare variable for distance_matrix
        VVD distance_matrix(N, VD(N, 0));

        // Read the distance_matrix from the input file into the variable
        general::ReadDoubleMatrix(inputfolder_parent + inputfolder_relative + inputfilename, distance_matrix, 0, N, 0, N);

        // Compute the corresponding normalized distance matrix
        NormalizeDistanceMatrix(distance_matrix, norm_distance_matrix);

        // Send the normalized distance matrix to each other rank
        for(int w = 1; w < world_size; w++)
        {
            MPI_Send(norm_distance_matrix.data(), N*N, MPI_DOUBLE, w, 0, MPI_COMM_WORLD);
        }
    }
    else
    {
        // Receive normalized distance matrix from rank 0
        MPI_Recv(norm_distance_matrix.data(), N*N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // Get all jobs i, i.e. compute the recurrence rate w.r.t. eps(i) = i * deps
    VI jobs_vector = CreateJobsVector(world_rank, world_size);
    size_t njobs_my_rank = jobs_vector.size();
    // For each rank and its jobs, save the result of the computation into field_distance_results
    VD recurrence_rate_results(njobs_my_rank, 0);

    // Compute the all recurrence rates based on jobs_vector
    double eps;
    double rr;
    for(size_t i = 0; i != njobs_my_rank; i++)
    {
        eps = jobs_vector[i] * deps;
        CalculateRecurrenceMatrix(norm_distance_matrix, eps, recurrence_matrix);
        rr = recurrence_rate(recurrence_matrix);
        recurrence_rate_results[i] = rr;
    }

    // Write the intermediate distance results 
    IO::WriteVectorResult(outputfolder_parent + outputfolder_relative + "recurrence_rate_" + std::to_string(world_rank) + ".txt", field_distance_results);

    // Rank 0 needs to compose the field_distance_matrix
    //if(world_rank == 0)
    //{
        // Save all distances from each rank into the following vector
        //VD total_field_distance_results(njobs, 0);
        // Calculate how the full problem (computing the distance_matrix) was split up into jobs
        //int njobs_per_rank_min = njobs / world_size;
        //int njobs_per_rank_max = njobs % world_size == 0 ? njobs_per_rank_min : njobs_per_rank_min + 1;
        //int min_rank = njobs % world_size;

        // Add own intermediate results to total resul vector
        //for(int i = 0; i < njobs_per_rank_max; i++)
        //{
        //    total_field_distance_results[i] = field_distance_results[i];
        //}

        // Receive all intermediate results, with vector lengths recv_result_size.
        // recv_offset is a variable which shifts the beginning of the receive buffer.
        //int recv_result_size;
        //int recv_offset = njobs_per_rank_max;
        //for(int w = 1; w < world_size; w++)
        //{
        //    recv_result_size = w < min_rank ? njobs_per_rank_max : njobs_per_rank_min;
        //    MPI_Recv(total_field_distance_results.data() + recv_offset, recv_result_size, MPI_DOUBLE, w, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //    recv_offset += recv_result_size;
        //}

        // Set the final field_distance_matrix from received results
        //VVD field_distance_matrix(nsteps, VD(nsteps, 0));
        //for(int i1 = 0, i = 0; i1 < nsteps - 1; i1++)
        //{
        //    for(int i2 = i1 + 1; i2 < nsteps; i2++, i++)
        //    {
        //        field_distance_matrix[i1][i2] = field_distance_matrix[i2][i1] = total_field_distance_results[i];
        //    }
        //}

        // Write the computed distance matrix based on fields into this filepathname
        //IO::WriteMatrixResult(outputfolder_parent + outputfolder_relative + outputfilename, field_distance_matrix);
    //}
    // Each other rank sends their intermediate results to rank 0
    //else
    //{
    //    MPI_Send(field_distance_results.data(), field_distance_results.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    //}

    // Finalize the MPI environment.
    MPI_Finalize();

    return 0;
}

// Calculate the recurrence rate rr for a give recurrence matrix
double recurrence_rate(VD &recurrence_matrix)
{
    double rr = 0;
    double sum = std::accumulate(recurrence_matrix.begin(), recurrence_matrix.end(), 0.0);
    int dim = N * N;
    rr = sum / dim;
    return rr;
}

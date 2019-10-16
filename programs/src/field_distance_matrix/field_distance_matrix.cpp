// Author:     Thomas Miethlinger B.Sc.
// Date:       16.10.2019
// Email:      t.miethlinger@protonmail.com
// University: Johannes Kepler University Linz, Austria

// STL
#include <iostream>
#include <string>
#include <cmath> // sqrt
#include <algorithm> // std::fill
#include <array>
#include <vector>
//#include <cctype>
//#include <cstdlib>
//#include <unistd.h>

// MPI
#include <mpi.h>

// Boost
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

// Internal classes
#include "../../includes/io_util.hpp"
#include "../../includes/general_util.hpp"

using std::size_t;
using std::string;
using std::array;
using std::vector;
using namespace boost::program_options;

typedef vector<int> VI;
typedef vector<double> VD;
typedef vector<VD> VVD;
typedef array<int, 2> AI2;

// General constants
const double precision = 0.000000000001;
const int digits = 3;

// Variables for the processing of the input dataset
// They are set by the programm arguments argv
int N; // number of particles
int imin;
int imax;
int istep;
int nsteps;
int njobs;

// Geometry specification
double xmin = -0.04;
double xmax = 0.04;
size_t nx;
double dx;

double ymin = -0.0075;
double ymax = 0.0075;
size_t ny;
double dy; // 0.00375

double zmin = 0.0;
double zmax = 0.25;
size_t nz;
double dz; // 0.00625

// Calculation
int get_index_from_coordinate(double x, double min, double d);
void calc_count_grid(VI& countgrid, VVD& R);
double L2_metric(VI& count_grid_a, VI& count_grid_b);

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
    ("N", value<int>()->default_value(1), "Particle number N")
    ("imin", value<int>()->default_value(0), "First time index")
    ("imax", value<int>()->default_value(1), "Last time index")
    ("istep", value<int>()->default_value(1), "Time index step size")
    ("nx", value<size_t>()->default_value(1), "Number of bins in x-direction")
    ("ny", value<size_t>()->default_value(1), "Number of bins in y-direction")
    ("nz", value<size_t>()->default_value(1), "Number of bins in z-direction")
    ("inputfolder_relative", value<string>()->default_value("Liggghts/"), "Relative input foldername")
    ("inputfilenamepart", value<string>()->default_value("dump"), "Input file name part")
    ("outputfolder_relative", value<string>()->default_value("Liggghts/"), "Relative output foldername");
    variables_map vm;
    store(parse_command_line(argc, argv, desc_commandline), vm);
    N = vm["N"].as<int>();
    imin = vm["imin"].as<int>();
    imax = vm["imax"].as<int>();
    istep = vm["istep"].as<int>();
    nsteps = (imax - imin) / istep + 1;
    njobs = (nsteps*nsteps - nsteps)/2;
    nx = vm["nx"].as<size_t>();
    dx = (xmax - xmin) / nx;
    ny = vm["ny"].as<size_t>();
    dy = (ymax - ymin) / ny;
    nz = vm["nz"].as<size_t>();
    dz = (zmax - zmin) / nz;
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
    // Define a vector which transforms from i to {i_1, i_2}
    vector<AI2> total_jobs_vector(njobs);
    for(int i1 = 0, i = 0; i1 < nsteps - 1; i1++)
    {
        for(int i2 = i1 + 1; i2 < nsteps; i2++, i++)
        {
            total_jobs_vector[i] = {i1, i2};
        }
    }
    // For each rank and its jobs, save the result of the computation into field_distance_results
    VD field_distance_results(njobs_my_rank, 0);

    // Declare two vectors for two states a, b, which count the number of particles per bin
    VI count_grid_a(nx * ny * nz, 0);
    VI count_grid_b(nx * ny * nz, 0);

    // Initialize particle vectors R_a and R_b
    VVD Ra(N, VD(3, 0));
    VVD Rb(N, VD(3, 0));

    // Variables for the indices of each job
    int i1, i2, i1_prev, i2_prev;
    i1 = i2 = i1_prev = i2_prev = -1;
    std::array<int, 2> job_index_arr;
    for(size_t i = 0; i != njobs_my_rank; i++)
    {
        job_index_arr = total_jobs_vector[jobs_vector_my_rank[i]];

        // Files may be saved per time-step with an offset imin and scaling istep.
        i1 = job_index_arr[0] * istep + imin;
        // In order to avoid two-times reading of the same file
        if(i1 != i1_prev)
        {
            // Reset count_grid_a
            std::fill(count_grid_a.begin(), count_grid_a.end(), 0);
            // Read input dump file
            if (!io_util::read_positions_3D(inputfolder_parent + inputfolder_relative + inputfilenamepart + std::to_string(i1), N, Ra))
                return 1;
            // Bin particles
            calc_count_grid(count_grid_a, Ra);
        }

        i2 = job_index_arr[1] * istep + imin;
        if(i2 != i2_prev)
        {
            std::fill(count_grid_b.begin(), count_grid_b.end(), 0);
            if (!io_util::read_positions_3D(inputfolder_parent + inputfolder_relative + inputfilenamepart + std::to_string(i2), N, Rb))
                return 1;
            calc_count_grid(count_grid_b, Rb);
        }

        // Compute L_2 metric
        field_distance_results[i] = L2_metric(count_grid_a, count_grid_b);
        i1_prev = i1;
        i2_prev = i2;
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

// Compute the bin index of coordinate x based on min. coordinate x and bin width d
int get_index_from_coordinate(double x, double min, double d)
{
    return static_cast<int>(std::floor((x - min)/d));
}

void calc_count_grid(VI& countgrid, VVD& R)
{
    int ix, iy, iz;
    for(int i = 0; i < N; i++)
    {
        ix = get_index_from_coordinate(R[i][0], xmin, dx);
        iy = get_index_from_coordinate(R[i][1], ymin, dy);
        iz = get_index_from_coordinate(R[i][2], zmin, dz);
        countgrid[(iz * ny + iy) * nx + ix]++;
    }
}

double L2_metric(VI& count_grid_a, VI& count_grid_b)
{
    double sum = 0.0;

    double d = 0.0;
    for(size_t iz = 0, i = 0; iz < nz; iz++)
    {
        for(size_t iy = 0; iy < ny; iy++)
        {
            for(size_t ix = 0; ix < nx; ix++, i++)
            {
                d = count_grid_a[i] - count_grid_b[i];
                sum += d*d;
            }
        }
    }

    return sqrt(sum);
}

// Determine the jobs per rank
/*vector<AI2> CreateJobsVector(int world_rank, int world_size)
{
    vector<AI2> job_vector;
    int njobs_per_rank_min = njobs / world_size;
    int njobs_per_rank_max = njobs % world_size == 0 ? njobs_per_rank_min : njobs_per_rank_min + 1;

    int min_rank = njobs % world_size;
    int ijobmin, ijobmax;
    if(world_rank < min_rank)
    {
        job_vector.resize(njobs_per_rank_max);
        ijobmin = world_rank * njobs_per_rank_max;
        ijobmax = ijobmin + njobs_per_rank_max - 1;
    }
    else
    {
        job_vector.resize(njobs_per_rank_min);
        ijobmin = min_rank * njobs_per_rank_max + (world_rank - min_rank) * njobs_per_rank_min;
        ijobmax = ijobmin + njobs_per_rank_min - 1;
    }

    for(int i1 = 0, i = 0, j = 0; i1 < nsteps - 1; i1++)
    {
        for(int i2 = i1 + 1; i2 < nsteps; i2++, i++)
        {
            if(ijobmin <= i && i <= ijobmax)
            {
                job_vector[j] = {i1, i2};
                j++;
            }
        }
    }

    return job_vector;
}*/

// MPI
#include <mpi.h>

// STL
#include <iostream>
#include <string>

// Math, Algorithm and Datastructures
#include <cmath>
#include <array>
#include <vector>

// Argument parsing
#include <cctype>
#include <cstdlib>
#include <unistd.h>

// Boost
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

// Internal classes
#include "ParticleIO.hpp"

typedef std::vector<double> VD;
typedef std::vector<VD> VVD;
typedef std::array<int, 2> AI2;
using std::size_t;
using namespace boost::program_options

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
int nx;
double dx; // 0.01

double ymin = -0.0075;
double ymax = 0.0075;
int ny;
double dy; // 0.00375

double zmin = 0.0;
double zmax = 0.25;
int nz;
double dz; // 0.00625

// Calculation
int GetTotalBinNumber(double min, double max, double d);
int GetIndexFromCoordinate(double x, double min, double d);
void FillCountGrid(VD& countgrid, VVD& R);
double L2Metric(VD& count_grid_a, VD& count_grid_b);
void ResetCountGrid(VD& countgrid);
void PrintCountGrid(VD& countgrid);

// Jobs
std::vector<AI2> CreateJobsVector(int world_rank, int world_size);

int main(int argc, char * argv[])
{
    MPI_Init(&argc, &argv);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    options_description desc_commandline;
    desc_commandline.add_options()
    ("N", value<int>()->default_value(1), "Particle number N")
    ("imin", value<int>()->default_value(0), "First time index")
    ("imax", value<int>()->default_value(1), "Last time index")
    ("istep", value<int>()->default_value(1), "Time index step size")
    ("nx", value<int>()->default_value(1), "Number of bins in x-direction")
    ("ny", value<int>()->default_value(1), "Number of bins in y-direction")
    ("nz", value<int>()->default_value(1), "Number of bins in z-direction");
    variables_map vm;
    store(parse_command_line(argc, argv, desc_commandline), vm);

    N = vm["N"].as<int>();
    imin = vm["imin"].as<int>();
    imax = vm["imax"].as<int>();
    istep = vm["istep"].as<int>();
    nsteps = (imax - imin + 1) / istep;
    njobs = (n_steps*n_steps - n_steps)/2;
    nx = vm["nx"].as<int>();
    dx = (xmax - xmin) / nx;
    ny = vm["ny"].as<int>();
    dy = (ymax - ymin) / ny;
    nz = vm["nz"].as<int>();
    dz = (zmax - zmin) / nz;

    // if(world_rank == 0)
    //     particle_io::CreateWorkingDirectory();

    std::vector<AI2> jobs_vector = CreateJobsVector(world_rank, world_size);

    VVD recurrenceresults(njobs, VD(3, 0));

    nx = GetTotalBinNumber(xmin, xmax, dx);
    ny = GetTotalBinNumber(ymin, ymax, dy);
    nz = GetTotalBinNumber(zmin, zmax, dz);

    VD count_grid_a(nx * ny * nz, 0):
    VD count_grid_b(nx * ny * nz, 0);

    // Initialize variables
    VVD Ra(N, VD(3, 0));
    VVD Rb(N, VD(3, 0));

    int i1 = i2 = i1_prev = i2_prev = -1;
    int t1, t2;
    std::array<int, 2> indices = {-1, -1};
    for(std::vector<AI2>::iterator it = jobs_vector.begin(); it != jobs_vector.end(); it++)
    {
        indices = *it;

        i1 = indices[0];
        recurrenceresults[ctr][0] = i1;
        if(i1 != i1_prev)
        {
            ResetCountGrid(count_grid_a);
            t1 = i1 * istep;
            if (!io.ReadLiggghtsDumpFull3D(Ra, t1))
                return 1;
            FillCountGrid(count_grid_a, Ra);
        }

        i2 = indices[1];
        recurrenceresults[ctr][1] = i2;
        if(i2 != i2_prev)
        {
            ResetCountGrid(count_grid_b);
            t2 = i2 * istep;
            if (!io.ReadLiggghtsDumpFull3D(Rb, t2))
                return 1;
            FillCountGrid(count_grid_b, Rb);
        }

        recurrenceresults[ctr][2] = L2Metric(count_grid_a, count_grid_b);
        i1_prev = i1;
        i2_prev = i2;
    }

    io.WriteRecurrenceMatrixResult(recurrenceresults, world_rank);

    // Finalize the MPI environment.
    MPI_Finalize();

    return 0;
}

int GetTotalBinNumber(double min, double max, double d)
{
    return static_cast<int>((max - min)/d + precision);
}

int GetIndexFromCoordinate(double x, double min, double d)
{
    return static_cast<int>(std::floor((x - min)/d));
}

void FillCountGrid(VD& countgrid, VVD& R)
{
    int ix, iy, iz;
    for(int i = 0; i < N; i++)
    {
        ix = GetIndexFromCoordinate(R[i][0], xmin, dx);
        iy = GetIndexFromCoordinate(R[i][1], ymin, dy);
        iz = GetIndexFromCoordinate(R[i][2], zmin, dz);
        countgrid[(iz * ny + iy) * nx + ix]++;
    }
}

double L2Metric(VD& count_grid_a, VD& count_grid_b)
{
    double sum = 0.0;

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

void ResetCountGrid(VD& countgrid)
{
    for(size_t iz = 0, i = 0; iz < nz; iz++)
    {
        for(size_t iy = 0; iy < ny; iy++)
        {
            for(size_t ix = 0; ix < nx; ix++, i++)
            {
                countgrid[i] = 0.0;
            }
        }
    }
}

void PrintCountGrid(VD& countgrid)
{
    for(size_t iz = 0, i = 0; iz < nz; iz++)
    {
        for(size_t iy = 0; iy < ny; iy++)
        {
            for(size_t ix = 0; ix < nx; ix++, i++)
            {
                std::cout << countgrid[i] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}

// Jobs
std::vector<AI2> CreateJobList(int world_rank, int world_size)
{
    std::vector<AI2> job_vector;
    int njobs_per_rank_min = njobs / world_size;
    int njobs_per_rank_max = njobs % world_size == 0 ? jobs_per_rank_min : jobs_per_rank_min + 1;

    int min_rank = njobs % world_size;
    if(world_rank < min_rank)
    {
        job_vector.resize(njobs_per_rank_max);
        ijobmin = world_rank * njobs_per_rank_max;
        ijobmax = ijobmax + njobs_per_rank_max - 1;
    }
    else
    {
        job_vector.resize(njobs_per_rank_min);
        ijobmin = min_rank * njobs_per_rank_max;
        ijobmax = ijobmax + njobs_per_rank_min - 1;
    }

    for(int i1 = 0, i = 0, j = 0; i1 < nsteps; i1++)
    {
        for(int i2 = i1 + 1; i2 <= nsteps; i2++, i++)
        {
            if(ijobmin <= i && i <= ijobmax)
            {
                job_vector[j] = {i1, i2};
                j++;
            }
        }
    }

    return job_vector;
}

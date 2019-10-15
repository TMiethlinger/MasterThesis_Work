// FIRST PARALLEL VERSION!

// MPI
#include <mpi.h>

// STL
#include <iostream>
#include <string>

// Math, Algorithm and Datastructures
#include <cmath>
#include <vector>

// Argument parsing
#include <cctype>
#include <cstdlib>
#include <unistd.h>

// Boost
#include <boost/algorithm/string.hpp>

// Internal classes
#include "ParticleIO.hpp"

// General constants
const double precision = 0.000000001;
const int digits = 3;
const int dim = 3;

// Constants and variables for this input dataset
int dtstep;

// Variables for the processing of the input dataset
// They are set by the programm arguments argv
int N; // number of particles
int istep_min;
int istep_max;

// Geometry specification
double xmin = -0.04;
double xmax = 0.04;
double dx = 0.01;
double ymin = -0.0075;
double ymax = 0.0075;
double dy = 0.00375;
double zmin = 0.0;
double zmax = 0.25;
double dz = 0.00625;

// Calculation
int GetTotalBinNumber(double min, double max, double d);
int GetIndexFromCoordinate(double x, double min, double d);
void FillCountGrid(std::vector<std::vector<std::vector<double>>>& countgrid, std::vector<std::vector<double>>& R);
double L2Metric(std::vector<std::vector<std::vector<double>>>& countgridi, std::vector<std::vector<std::vector<double>>>& countgridj);
void ResetCountGrid(std::vector<std::vector<std::vector<double>>>& countgrid);
void PrintCountGrid(std::vector<std::vector<std::vector<double>>>& countgrid);

// Jobs
std::list<std::pair<int,int>> CreateJobList(int world_rank, int world_size);

int main(int argc, char * argv[])
{
	int world_size;
	int world_rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	char c;
	while ((c = getopt(argc, argv, "N:i:j:d:")) != -1)
	{
		switch (c)
		{
			case 'N':
				N = std::stoi(optarg);
				break;
			case 'i':
				istep_min = std::stoi(optarg);
				break;
			case 'j':
				istep_max = std::stoi(optarg);
				break;
			case 'd':
				dtstep = std::stoi(optarg);
				break;
			default:
				abort();
		}
	}

	ParticleIO io(N, istep_min, istep_max, dtstep);
	io.CreateWorkingDirectory(world_rank);

	std::list<std::pair<int,int>> joblist = CreateJobList(world_rank, world_size);

	std::vector<std::vector<double>> recurrenceresults(joblist.size(), std::vector<double>(3, 0.0));

	int nbinsx = GetTotalBinNumber(xmin, xmax, dx);
	int nbinsy = GetTotalBinNumber(ymin, ymax, dy);
	int nbinsz = GetTotalBinNumber(zmin, zmax, dz);

	std::vector<std::vector<std::vector<double>>> countgridi(nbinsx, std::vector<std::vector<double>>(nbinsy, std::vector<double>(nbinsz, 0)));
	std::vector<std::vector<std::vector<double>>> countgridj(nbinsx, std::vector<std::vector<double>>(nbinsy, std::vector<double>(nbinsz, 0)));

	// Initialize variables
	std::vector<std::vector<double>> Ri(N, std::vector<double>(dim));
	std::vector<std::vector<double>> Rj(N, std::vector<double>(dim));

	int i, j, ctr;
	ctr = 0;
	int ti, tj;
	for(std::list<std::pair<int,int>>::iterator it = joblist.begin(); it != joblist.end(); it++, ctr++)
    {
		if(ctr == 0)
		{
			i = (*it).first;
			recurrenceresults[ctr][0] = i;
			ti = i * dtstep;
			if (!io.ReadLiggghtsDumpFull3D(Ri, ti))
			{
				return 1;
			}
			FillCountGrid(countgridi, Ri);

			j = (*it).second;
			recurrenceresults[ctr][1] = j;
	        tj = j * dtstep;
			if (!io.ReadLiggghtsDumpFull3D(Rj, tj))
			{
				return 1;
			}
			FillCountGrid(countgridj, Rj);

			recurrenceresults[ctr][2] = L2Metric(countgridi, countgridj);
		}
		else
		{
			i = (*it).first;
			recurrenceresults[ctr][0] = i;
			if(i != recurrenceresults[ctr - 1][0])
			{
				ResetCountGrid(countgridi);
				ti = i * dtstep;
				if (!io.ReadLiggghtsDumpFull3D(Ri, ti))
				{
					return 1;
				}
				FillCountGrid(countgridi, Ri);
			}

			j = (*it).second;
			recurrenceresults[ctr][1] = j;
			if(j != recurrenceresults[ctr - 1][1])
			{
				ResetCountGrid(countgridj);
				tj = j * dtstep;
				if (!io.ReadLiggghtsDumpFull3D(Rj, tj))
				{
					return 1;
				}
				FillCountGrid(countgridj, Rj);
			}

			recurrenceresults[ctr][2] = L2Metric(countgridi, countgridj);
		}
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

void FillCountGrid(std::vector<std::vector<std::vector<double>>>& countgrid, std::vector<std::vector<double>>& R)
{
	int nx, ny, nz;
	for(auto v1 : R)
	{
		nx = GetIndexFromCoordinate(v1[0], xmin, dx);
		ny = GetIndexFromCoordinate(v1[1], ymin, dy);
		nz = GetIndexFromCoordinate(v1[2], zmin, dz);
		countgrid[nx][ny][nz]++;
	}
}

double L2Metric(std::vector<std::vector<std::vector<double>>>& countgridi, std::vector<std::vector<std::vector<double>>>& countgridj)
{
	double sum = 0.0;
	double d = 0.0;
	int n1 = countgridi.size();
	int n2 = countgridi[0].size();
	int n3 = countgridi[0][0].size();
	for(int i = 0; i < n1; i++)
	{
		for(int j = 0; j < n2; j++)
		{
			for(int k = 0; k < n3; k++)
			{
				d = countgridi[i][j][k] - countgridj[i][j][k];
				sum += d*d;
			}
		}
	}

	return sqrt(sum);
}

void ResetCountGrid(std::vector<std::vector<std::vector<double>>>& countgrid)
{
	int n1 = countgrid.size();
	int n2 = countgrid[0].size();
	int n3 = countgrid[0][0].size();
	for(int i = 0; i < n1; i++)
	{
		for(int j = 0; j < n2; j++)
		{
			for(int k = 0; k < n3; k++)
			{
				countgrid[i][j][k] = 0.0;
			}
		}
	}
}

void PrintCountGrid(std::vector<std::vector<std::vector<double>>>& countgrid)
{
	for(auto v2 : countgrid)
	{
		for(auto v1 : v2)
		{
			for(double d : v1)
			{
				std::cout << d << " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
}

// Jobs
std::list<std::pair<int,int>> CreateJobList(int world_rank, int world_size)
{
    // Return object: list of pairs of (i,j). i and j are indices for the timesteps which are to be compared.
    std::list<std::pair<int,int>> joblist;

    // Count job list
    int n_jobs = 0;
    for(int i = istep_min; i < istep_max; i++)
    {
        for(int j = i + 1; j <= istep_max; j++)
        {
            n_jobs++;
        }
    }

	// If the number of jobs is equal to the world_size:
	// -> joblist.size() == 1
    if(n_jobs == world_size)
    {
		int k = 0;
       	for(int i = istep_min; i < istep_max; i++)
        {
            for(int j = i + 1; j <= istep_max; j++, k++)
            {
                if(k == world_rank)
                {
					joblist.push_back(std::make_pair(i, j));
				}
			}
		}
	}
	// If the number of jobs is smaller than the world_size:
	// joblist.size() = 0 or 1
    else if(n_jobs < world_size)
    {
		if(world_rank < n_jobs) // e.g. {1, ..., 7} < 8 -> do something.
		{
			int k = 0;
	       	for(int i = istep_min; i < istep_max; i++)
	        {
	            for(int j = i + 1; j <= istep_max; j++, k++)
	            {
	                if(k == world_rank)
	                {
						joblist.push_back(std::make_pair(i, j));
					}
				}
			}
		}
		else // e.g. {8, 9, ... } >= 8 -> do nothing
		{
		}
	}
	// If the number of jobs is larger than the world_size:
	// joblist.size() >= 1
    else if(world_size < n_jobs)
    {
	    // If the jobs can be evenly distributed:
		// joblist.size() = n_jobs / world_size
	    if(n_jobs % world_size == 0)
	    {
	        int jobs = n_jobs / world_size;
	        int k = 0;
	        int k_min = jobs * world_rank;
	        int k_max = (jobs * (world_rank + 1)) - 1;
	        for(int i = istep_min; i < istep_max; i++)
	        {
	            for(int j = i + 1; j <= istep_max; j++, k++)
	            {
	                if(k_min <= k && k <= k_max)
	                {
	                    joblist.push_back(std::make_pair(i, j));
	                }
	            }
	        }
	    }
		// Else: non-trivial number of jobs for each node.
		else
		{
			// e.g. n_jobs = 43 jobs on world_size = 5 processes
    		// perfect matching: 9,9,9,8,8 -> 0...8, 9...17, 18...26, 27...34, 35...42
			// Algorithm: jobs_min := 43 / 5 = 8; jobs_max := jobs_min + 1 = 9.
			// 1.) n_jobs - world_size * (n_jobs / world_size) = 43 - 40 = 3 =: border -> additional jobs on rank0, rank1 and rank2
			// 2.) IF(world_rank < border): k_min = world_rank * jobs_max & k_max = (world_rank + 1) * jobs_max - 1. Example: world_rank = 2. k_min = 18 & k_max = 26.
			// 3.) ELSE: k_min = border * jobs_max + (world_rank - border) * jobs_min & k_max = border * jobs_max + ((world_rank + 1 - border) * jobs_min - 1)
			// Example a.): world_rank = 3. k_min = 3 * 9 + (3 - 3) * 8 = 27 & k_max = 3 * 9 + (1 * 8 - 1) = 34
			// Example b.): world_rank = 4. k_min = 3 * 9 + 1 * 8 = 35 & k_max = 3 * 9 + 2 * 8 - 1 = 27 + 16 - 1 = 42

			int jobs_min = n_jobs / world_size;
			int jobs_max = jobs_min + 1;
			int border = n_jobs - world_size * (n_jobs / world_size);
			int k_min;
			int k_max;
			if(world_rank < border)
			{
				k_min = world_rank * jobs_max;
				k_max = (world_rank + 1) * jobs_max - 1;
			}
			else
			{
				k_min = border * jobs_max + (world_rank - border) * jobs_min;
				k_max = border * jobs_max + ((world_rank + 1 - border) * jobs_min - 1);
			}

	        int k = 0;
	        for(int i = istep_min; i < istep_max; i++)
	        {
	            for(int j = i + 1; j <= istep_max; j++, k++)
	            {
	                if(k_min <= k && k <= k_max)
	                {
	                    joblist.push_back(std::make_pair(i, j));
	                }
	            }
	        }
		}
	}

	return joblist;
}

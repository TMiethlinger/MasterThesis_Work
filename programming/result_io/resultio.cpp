/// HOW TO USE:
/// argv: ./resultio datasetfolderpath relinputfolderpath reloutputfolderpath outputfilename imin imax di
/// ./resultio /home/k3501/k354524/MasterThesis/Programming/HungarianMethod/cpp/Production/8_0_lammps_2d_pbc_G/results/100_0_0_11.180340_1_100000_500000_10/ elements/ analysis/ distancematrix.txt 100000 500000 10
/// HOW TO COMPILE:
/// g++ -std=c++17 -w -O3 -o resultio resultio.cpp -lstdc++fs

/* STL */
// IO and files
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sys/stat.h>
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem::v1;

// String
#include <string>
// Containers
#include <vector>

/// External libraries
// Boost
#include <boost/algorithm/string.hpp>

// general
std::string datasetfolderpath;
// input
std::string relinputfolderpath;
std::string absinputfolderpath;
// output
std::string reloutputfolderpath;
std::string absoutputfolderpath;
std::string outputfilename;
std::string outputfilepathname;

int imin;
int imax;
int di;
int n;

// IO
bool ReadCostFile(std::string inputfilepathname, std::vector<std::vector<double>>& costmatrix);
bool WriteCostMatrix(std::string outputfilepathname, std::vector<std::vector<double>> &costmatrix);

int main(int argc, char * argv[])
{
    // Set variables from arguments list
    datasetfolderpath = argv[1];
    relinputfolderpath = argv[2];
    reloutputfolderpath = argv[3];
    outputfilename = argv[4];
    imin = std::stoi(argv[5]);
    imax = std::stoi(argv[6]);
    di = std::stoi(argv[7]);

    absinputfolderpath = datasetfolderpath + relinputfolderpath;
    absoutputfolderpath = datasetfolderpath + reloutputfolderpath;
    mkdir(absoutputfolderpath.c_str(), 0777);
    outputfilepathname = absoutputfolderpath + outputfilename;
    n = (imax - imin) / di + 1;

    std::vector<std::vector<double>> costmatrix(n, std::vector<double>(n, 0.0));

    for (const auto & entry : fs::directory_iterator(absinputfolderpath))
    {
        std::cout << entry.path() << std::endl;
        if(!ReadCostFile(entry.path(), costmatrix))
        {
            std::cout << "!ReadCostFile(" << entry.path() << ", costmatrix)" << std::endl;
            return 1;
        }
    }

    if(!WriteCostMatrix(outputfilepathname, costmatrix))
    {
        std::cout << "!WriteCostMatrix(" << outputfilepathname << ", costmatrix)" << std::endl;
        return 1;
    }

    return 0;
}

// IO
bool ReadCostFile(std::string inputfilepathname, std::vector<std::vector<double>>& costmatrix)
{
    int i, j;
    double cost;

    std::ifstream filestream;
    filestream.open(inputfilepathname);
    if (filestream.is_open())
    {
        std::string line;
        std::vector<std::string> parts;
        for (int linectr = 0; std::getline(filestream, line); linectr++)
        {
            boost::split(parts, line, boost::is_any_of(" "));
	    
            i = (std::stoi(parts[0]) - imin) / di;
            j = (std::stoi(parts[1]) - imin) / di;
            cost = std::stod(parts[2]);

	    if(i < n && j < n)
	    {
	        costmatrix[i][j] = cost;
		costmatrix[j][i] = costmatrix[i][j];
	    }
        }
        filestream.close();
        return true;
    }
    else
        return false;
}

bool WriteCostMatrix(std::string outputfilepathname, std::vector<std::vector<double>> &costmatrix)
{
    std::ofstream filestream;
    filestream.open(outputfilepathname);
    if (filestream.is_open())
    {
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
            {
                filestream << costmatrix[i][j];
                if(j < n - 1)
                    filestream << " ";
            }
            if(i < n - 1)
                filestream << std::endl;
        }
        filestream.close();
        return true;
    }
    else
        return false;
}

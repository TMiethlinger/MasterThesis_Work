// Author: Thomas Miethlinger BSc.
// thomas.miethlinger@gmail.com

// STL
#include <iostream>
#include <iomanip>
#include <fstream>

#include <chrono>
#include <thread>
#include <sys/stat.h>

#include <string>

#include <vector>

// Boost
#include <boost/algorithm/string.hpp>

#pragma once

typedef std::vector<double> VD;
typedef std::vector<VD> VVD;

namespace particle_io
{
    void CreateWorkingDirectory(char *folderpath)
    {
        mkdir(folderpath, 0777);
    }

    bool ReadPositions3D(std::string filepathname, int N, VVD &R)
    {
        bool result = true;

        std::ifstream filestream;
        filestream.open(filepathname);
        if (filestream.is_open())
        {
            std::string line;
            std::vector<std::string> parts;
            for (int i = -9; std::getline(filestream, line) && i < N; i++)
            {
                if (i >= 0)
                {
                    boost::split(parts, line, boost::is_any_of(" "));
                    VD r(3);
                    for (int j = 0; j < 3; j++)
                    {
                        r[j] = std::stod(parts[j + 1]);
                    }
                    R[i] = vec;
                }
            }
            filestream.close();
        }
        else
        {
            std::cout << "Error: Could not read file " << inputfilepathname << "!" << std::endl;
            result = false;
        }

        return result;
    }

    bool WriteRecurrenceMatrixResult(std::string filepathname, VVD& recurrence_matrix)
    {
        bool success = false;

        std::ofstream filestream;
        int ctrfail = 0;
        do
        {
            filestream.open(filepathname);
            if(filestream.is_open())
            {
                for(int i = 0; i < recurrence_matrix.size(); i++)
                {
                    for(int j = 0; j < recurrence_matrix[i].size(); j++)
                    {
                        filestream << recurrence_matrix[i][j];
                        if(j < recurrence_matrix[i].size() - 1)
                            filestream << " ";
                    }
                    filestream << std::endl;
                }

                success = true;
                filestream.close();
            }
            else
            {
                success = false;
                ctrfail++;
                std::this_thread::sleep_for(std::chrono::milliseconds(1000));
            }
        } while(!success && ctrfail < 10);
        if(!success)
        {
            std::cout << "Error: Could not write to file " << outputfilepathname << "!" << std::endl;
        }

        return success;
    }
};

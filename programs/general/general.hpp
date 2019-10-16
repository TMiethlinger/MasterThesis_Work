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

typedef std::vector<int> VI;
typedef std::vector<VI> VVI;
typedef std::vector<double> VD;
typedef std::vector<VD> VVD;
using std::size_t;
using std::string;

namespace general
{
    void CreateWorkingDirectory(const char *folderpath)
    {
        if (mkdir(folderpath, 0777) == -1 && errno != 17) 
            std::cerr << "Error: " << strerror(errno) << std::endl; 
    }

    bool ReadDoubleMatrix(std::string filepathname, VVD &m, string del, int row_min, int row_ctr, int col_min, int col_ctr)
    {
        bool result = true;

        m.resize(row_ctr, VD(col_ctr, 0));

        std::ifstream filestream;
        filestream.open(filepathname);
        if (filestream.is_open())
        {
            std::string line;
            std::vector<std::string> parts;
            for (int i = 0, ctr = 0; std::getline(filestream, line) && ctr < row_ctr; i++)
            {
                if (row_min <= i)
                {
                    boost::split(parts, line, boost::is_any_of(del));
                    VD v(col_ctr);
                    for (int j = 0; j < col_ctr; j++)
                    {
                        v[j] = std::stod(parts[j + col_min]);
                    }
                    m[ctr] = v;
                }
            }
            filestream.close();
        }
        else
        {
            std::cout << "Error: Could not read file " << filepathname << "!" << std::endl;
            result = false;
        }

        return result;
    }

    bool WriteVectorResult(std::string filepathname, VD& v)
    {
        bool success = false;

        std::ofstream filestream;
        int ctrfail = 0;
        do
        {
            filestream.open(filepathname);
            if(filestream.is_open())
            {
                for(size_t i = 0; i < v.size(); i++)
                {
                    filestream << v[i] << std::endl;
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
            std::cout << "Error: Could not write to file " << filepathname << "!" << std::endl;
        }

        return success;
    }

    bool WriteMatrixResult(std::string filepathname, VVD& m)
    {
        bool success = false;

        std::ofstream filestream;
        int ctrfail = 0;
        do
        {
            filestream.open(filepathname);
            if(filestream.is_open())
            {
                for(size_t i = 0; i < m.size(); i++)
                {
                    for(size_t j = 0; j < m[i].size(); j++)
                    {
                        filestream << m[i][j];
                        if(j < m[i].size() - 1)
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
            std::cout << "Error: Could not write to file " << filepathname << "!" << std::endl;
        }

        return success;
    }
};

// Author:     Thomas Miethlinger B.Sc.
// Date:       16.10.2019
// Email:      t.miethlinger@protonmail.com
// University: Johannes Kepler University Linz, Austria

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

namespace io_util
{
    // General IO functions

    // Create directory or write error unless the directory already exists
    void create_working_directory(const char *folderpath)
    {
        if (mkdir(folderpath, 0777) == -1 && errno != 17) // errno == 17 -> directory already exists
            std::cerr << "Error: " << strerror(errno) << std::endl; 
    }

    bool read_double_matrix(std::string filepathname, std::vector<std::vector<double>> &matrix, std::string del, int row_min, int row_ctr, int col_min, int col_ctr)
    {
        bool result = true;

        matrix.resize(row_ctr, std::vector<double>(col_ctr, 0));
        std::vector<double> row(col_ctr);

        std::ifstream filestream;
        filestream.open(filepathname);
        if (filestream.is_open())
        {
            std::string line;
            std::vector<std::string> parts;
            for (int i1 = 0, i2 = 0; std::getline(filestream, line) && i2 < row_ctr; i1++)
            {
                if (row_min <= i1)
                {
                    boost::split(parts, line, boost::is_any_of(del));
                    for (int j = 0; j < col_ctr; j++)
                    {
                        row[j] = std::stod(parts[j + col_min]);
                    }
                    matrix[i2] = row;
                    i2++;
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

    bool write_vector_result(std::string filepathname, std::vector<double>& vector)
    {
        bool success = false;

        std::ofstream filestream;
        filestream.open(filepathname);
        if(filestream.is_open())
        {
            for(std::size_t i = 0; i < vector.size(); i++)
            {
                filestream << vector[i] << std::endl;
            }

            success = true;
            filestream.close();
        }
        else
        {
            success = false;
            std::cout << "Error: Could not write to file " << filepathname << "!" << std::endl;
        }

        return success;
    }

    bool write_matrix_result(std::string filepathname, std::vector<std::vector<double>>& matrix)
    {
        bool success = false;

        std::ofstream filestream;
        filestream.open(filepathname);
        if(filestream.is_open())
        {
            for(std::size_t i = 0; i < matrix.size(); i++)
            {
                for(std::size_t j = 0; j < matrix[i].size(); j++)
                {
                    filestream << matrix[i][j];
                    if(j < matrix[i].size() - 1)
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
            std::cout << "Error: Could not write to file " << filepathname << "!" << std::endl;
        }

        return success;
    }



    // IO functions specifically to read/write particle data

    bool read_positions_3D(std::string filepathname, int N, std::vector<std::vector<double>> &R)
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
                    std::vector<double> r(3);
                    for (int j = 0; j < 3; j++)
                    {
                        r[j] = std::stod(parts[j + 1]);
                    }
                    R[i] = r;
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
};

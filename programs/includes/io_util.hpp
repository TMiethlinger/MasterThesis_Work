// Author:     Thomas Miethlinger B.Sc.
// Date:       16.10.2019
// Email:      t.miethlinger@protonmail.com
// University: Johannes Kepler University Linz, Austria

// STL
#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <list>
#include <thread>
#include <sys/stat.h>
#include <string>
#include <vector>
#include <utility> // std::pair

// Boost
#include <boost/algorithm/string.hpp>

#pragma once

namespace io_util
{
    // General IO functions

    // Create directory or write error unless the directory already exists
    void create_working_directory(const char *folderpath)
    {
        std::string folderpath_s(folderpath);
        if (mkdir(folderpath, 0777) == -1 && errno != 17) // errno == 17 -> directory already exists
            std::cerr << "Error: " << strerror(errno) << " (" << folderpath_s << ")" << std::endl; 
    }

    template <class T>
    void print_matrix_row_major_mode(std::vector<T> &matrix, int m, int n)
    {
        for(int i = 0; i < m; i++)
        {
            for(int j = 0; j < n; j++)
            {
                std::cout << matrix[i * n + j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    template <class T>
    void print_matrix(std::vector<std::vector<T>> &matrix)
    {
        for(std::size_t i = 0; i < matrix.size(); i++)
        {
            for(std::size_t j = 0; j < matrix[i].size(); j++)
            {
                std::cout << matrix[i][j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    template <class T>
    void print_vector(std::vector<T> &vec)
    {
        for(std::size_t i = 0; i < vec.size(); i++)
        {
            std::cout << vec[i] << std::endl;
        }
        std::cout << std::endl;
    }

    void print_adjlist(std::vector<std::vector<std::pair<int, double>>> &adjlist)
    {    
        for(size_t i = 0; i < adjlist.size(); i++)
        {
            for(size_t j = 0; j < adjlist[i].size(); j++)
            {
                std::cout << "(" << adjlist[i][j].first << ", " << adjlist[i][j].second << ") ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    bool read_double_matrix(std::string filepathname, std::vector<std::vector<double>> &matrix, std::string del, int row_min, int row_ctr, int col_min, int col_ctr)
    {
        bool result = true;

        if((std::size_t)row_ctr != matrix.size() || (std::size_t)col_ctr != matrix[0].size())
        {
            std::cout << "Error: Matrix dimensions do not fit! " << row_ctr << " " << matrix.size() << ", " << col_ctr << " " << matrix[0].size() << std::endl;
            return false;
        }

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

    bool read_double_matrix_row_major(std::string filepathname, std::vector<double> &matrix, std::string del, int row_min, int row_ctr, int col_min, int col_ctr)
    {
        bool result = true;

        if((std::size_t)(row_ctr * col_ctr) != matrix.size())
        {
            std::cout << "Error: Matrix dimensions do not fit! " << (row_ctr * col_ctr) << " " << matrix.size() << std::endl;
            return false;
        }

        std::ifstream filestream;
        filestream.open(filepathname);
        if (filestream.is_open())
        {
            std::string line;
            std::vector<std::string> parts;
            for (int i1 = 0, i2 = 0, i3 = 0; std::getline(filestream, line) && i2 < row_ctr; i1++)
            {
                if (row_min <= i1)
                {
                    boost::split(parts, line, boost::is_any_of(del));
                    for (int j = 0; j < col_ctr; j++, i3++)
                    {
                        matrix[i3] = std::stod(parts[j + col_min]);
                    }
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

    std::list<std::pair<int,int>> read_result_times_discrete_distance(std::string filepathname)
    {
        std::list<std::pair<int, int>> result_times;
        std::pair<int, int> row;
        std::ifstream filestream;
        filestream.open(filepathname);
        if(filestream.is_open())
        {
            std::string line;
            std::vector<std::string> parts;
            for(int i = 0; std::getline(filestream, line); i++)
            {
                boost::split(parts, line, boost::is_any_of(" "));
                row = std::make_pair(std::stoi(parts[0]), std::stoi(parts[1]));
                result_times.push_back(row);
            }
            filestream.close();
        }
        else
        {
            std::cout << "Error: Could not read file " << filepathname << "!" << std::endl;
        }

        return result_times;
    }

    template <class T>
    bool write_vector_result(std::string filepathname, std::vector<T>& vector)
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

    template <class T>
    bool write_matrix_result(std::string filepathname, std::vector<std::vector<T>>& matrix)
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

    template <class T>
    bool write_matrix_result_appending(std::string filepathname, std::vector<std::vector<T>>& matrix)
    {
        bool success = false;

        std::ofstream filestream;
        filestream.open(filepathname, std::ios_base::app);
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

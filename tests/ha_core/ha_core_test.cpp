// Author:     Thomas Miethlinger B.Sc.
// Date:       01.11.2019
// Email:      t.miethlinger@protonmail.com
// University: Johannes Kepler University Linz, Austria

// STL
#include <iostream>
#include <vector>
#include <random>
#include <functional>
#include <chrono> 

// Boost
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

// Internal classes
#include "../../programs/includes/io_util.hpp"
#include "../../programs/src/hungarian_algorithm/ha_core.hpp"
#include "../../programs/src/hungarian_algorithm/ha_cost.hpp"
#include "../../programs/src/hungarian_algorithm/ha_distance.hpp"

using std::size_t;
using std::string;
using std::vector;
using std::cout;
using std::endl;
using namespace std::chrono; 
using namespace boost::program_options;

typedef vector<int> VI;
typedef vector<double> VD;
typedef vector<VD> VVD;
typedef std::pair<int, double> PID;

int N;
int N_match;
int dim;
double l = 1.0;
double v = 1.0;

bool check_unique(VI &q_list);
int count_diff(int N, VI &q_mat, VI &q_list);
int count_unmatched(int N, VI &q_list);

int main(int argc, char * argv[])
{
    options_description desc_commandline;
    desc_commandline.add_options()
    ("N", value<int>()->default_value(10), "Number of vertices N")
    ("N_match", value<int>()->default_value(4), "Number of vertices N_match used for matching.")
    ("dim", value<int>()->default_value(6), "State vector dimension");
    variables_map vm;
    store(parse_command_line(argc, argv, desc_commandline), vm);
    N = vm["N"].as<int>();
    N_match = vm["N_match"].as<int>();
    dim = vm["dim"].as<int>();

    VVD Xa(N, VD(dim, 0.0));
    VVD Xb(N, VD(dim, 0.0));

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < dim; j++)
        {
            Xa[i][j] = dis(gen);
        }
    }
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < dim; j++)
        {
            Xb[i][j] = dis(gen);
        }
    }

    VVD cost_adjmatrix = ha_cost::create_costobject_adjmatrix(N, N_match, Xa, Xb, ha_distance::d_sum_3, l, v);
    vector<vector<PID>> cost_adjlist = ha_cost::create_costobject_adjlist_plain(N, N_match, Xa, Xb, ha_distance::d_sum_3, l, v);


    auto start_mat = high_resolution_clock::now();
    VI q_mat = ha_core::min_cost_matching_adjmatrix(N, cost_adjmatrix);
    auto duration_mat = duration_cast<nanoseconds>(high_resolution_clock::now() - start_mat);
    double c_mat = ha_cost::compute_total_cost_adjmatrix(N, cost_adjmatrix, q_mat);


    auto start_list_1 = high_resolution_clock::now();
    VI q_list_1 = ha_core::min_cost_matching_adjlist_1(N, cost_adjlist);
    auto duration_list_1 = duration_cast<nanoseconds>(high_resolution_clock::now() - start_list_1);
    double c_list_1 = ha_cost::compute_total_cost_adjlist(N, cost_adjlist, q_list_1);
    bool same_list_1 = c_mat == c_list_1 && q_mat == q_list_1;
    bool is_unique_list_1 = check_unique(q_list_1);
    int diff_list_1 = count_diff(N, q_mat, q_list_1);
    int unmatched_list_1 = count_unmatched(N, q_list_1);


    auto start_list_2 = high_resolution_clock::now();
    VI q_list_2 = ha_core::min_cost_matching_adjlist_3(N, cost_adjlist, N_match);
    auto duration_list_2 = duration_cast<nanoseconds>(high_resolution_clock::now() - start_list_2);
    double c_list_2 = ha_cost::compute_total_cost_adjlist(N, cost_adjlist, q_list_2);
    bool same_list_2 = c_mat == c_list_2 && q_mat == q_list_2;
    bool is_unique_list_2 = check_unique(q_list_2);
    int diff_list_2 = count_diff(N, q_mat, q_list_2);
    int unmatched_list_2 = count_unmatched(N, q_list_2);

    cout << "Number of left vertices N: " << N << endl;
    cout << "Number of right vertices per left vertex N_match: " << N_match << endl;

    //io_util::print_matrix(cost_adjmatrix);
    //io_util::print_adjlist(cost_adjlist);

    cout << "duration_mat [ms]: " << (duration_mat.count() / 1000000) << endl;
    cout << "duration_list_1 [ms]: " << (duration_list_1.count() / 1000000) << endl;
    cout << "duration_list_2 [ms]: " << (duration_list_2.count() / 1000000) << endl << endl;

    //io_util::print_vector(q_mat);
    //io_util::print_vector(q_list_1);
    //io_util::print_vector(q_list_2);

    cout << "cost_mat: " << c_mat << endl;
    cout << "c_list_1: " << c_list_1 << endl;
    cout << "c_list_2: " << c_list_2 << endl << endl;

    cout << "same_list_1: " << same_list_1 << endl;
    cout << "same_list_2: " << same_list_2 << endl << endl;

    cout << "is_unique_list_1: " << is_unique_list_1 << endl;
    cout << "is_unique_list_2: " << is_unique_list_2 << endl << endl;

    cout << "diff_list_1: " << diff_list_1 << endl;
    cout << "diff_list_2: " << diff_list_2 << endl << endl;

    cout << "unmatched_list_1: " << unmatched_list_1 << endl;
    cout << "unmatched_list_2: " << unmatched_list_2 << endl << endl;
}

bool check_unique(VI &q_list)
{
    auto it = std::unique(q_list.begin(), q_list.end());
    return (it == q_list.end());
}
int count_diff(int N, VI &q_mat, VI &q_list)
{
    int ctr = 0;
    for(int i = 0; i < N; i++)
    {
        if(q_mat[i] != q_list[i])
                ctr++;
    }
    return ctr;
}
int count_unmatched(int N, VI &q_list)
{
    int ctr = 0;
    for(int i = 0; i < N; i++)
    {
        if(-1 == q_list[i])
            ctr++;
    }
    return ctr;
}

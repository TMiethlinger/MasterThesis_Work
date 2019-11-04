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

for(int t = 0; t < 1; t++)
{

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
    //io_util::print_matrix(cost_adjmatrix);
    VI q_mat = ha_core::min_cost_matching_adjmatrix(N, cost_adjmatrix);
    //io_util::print_vector(q_mat);
    double c_mat = ha_cost::compute_total_cost_adjmatrix(N, cost_adjmatrix, q_mat);
    // cout << c << endl << endl << endl << endl;

    vector<vector<PID>> cost_adjlist = ha_cost::create_costobject_adjlist_plain(N, N_match, Xa, Xb, ha_distance::d_sum_3, l, v);
    //io_util::print_adjlist(cost_adjlist);
    auto start_2 = high_resolution_clock::now();
    VI q_list = ha_core::min_cost_matching_adjlist_1(N, cost_adjlist);
    auto stop_2 = high_resolution_clock::now();
    auto duration_2 = duration_cast<nanoseconds>(stop_2 - start_2);
    cout << "duration2:" << (duration_2.count()/1000000) << endl;
    //io_util::print_vector(q_list);
    double c_list = ha_cost::compute_total_cost_adjlist(N, cost_adjlist, q_list);
    // cout << c << endl << endl;
    bool same = c_mat == c_list && q_mat == q_list;
    if(true)
    {
        int count_diff = 0;
        int count_unmatched = 0;
        auto it = std::unique(q_list.begin(), q_list.end());
        bool wasUnique = (it == q_list.end());
        for(int i = 0; i < N; i++)
        { 
            if(q_mat[i] != q_list[i])
                count_diff++;

            if(-1 == q_list[i])
                count_unmatched++;
        }
        cout << "same: " << same << " count_unmatched:" << count_unmatched << " wasUnique:" << wasUnique << endl;
        cout << c_mat << " " << c_list << " " << (c_list - c_mat) << endl;
        cout << ha_cost::compute_total_cost_adjmatrix(N, cost_adjmatrix, q_list) << " " << ha_cost::compute_total_cost_adjlist(N, cost_adjlist, q_mat) << endl;
    }
}
}

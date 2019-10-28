// Author:     Thomas Miethlinger B.Sc.
// Date:       28.10.2019
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
#include "../../programs/src/hungarian_algorithm/ha_costobject.hpp"
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
int dim;
double l = 1.0;
double v = 1.0;

int main(int argc, char * argv[])
{
    options_description desc_commandline;
    desc_commandline.add_options()
    ("N", value<int>()->default_value(10), "Particle number N")
    ("dim", value<int>()->default_value(6), "State vector dimension");
    variables_map vm;
    store(parse_command_line(argc, argv, desc_commandline), vm);
    N = vm["N"].as<int>();
    dim = vm["N"].as<int>();

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

    VVD cost_adjmatrix = ha_costobject::create_costobject_adjmatrix(N, Xa, Xb, ha_distance::d_sum_3, l, v);
    io_util::print_matrix(cost_adjmatrix);

    int N_neighbours = N / 2;
    cost_adjmatrix = ha_costobject::create_costobject_adjmatrix(N, N_neighbours, Xa, Xb, ha_distance::d_sum_3, l, v);
    io_util::print_matrix(cost_adjmatrix);

    vector<vector<PID>> cost_adjlist = ha_costobject::create_costobject_adjlist_plain(N, N_neighbours, Xa, Xb, ha_distance::d_sum_3, l, v);

    for(size_t i = 0; i < cost_adjlist.size(); i++)
    {
        for(size_t j = 0; j < cost_adjlist[i].size(); j++)
        {
            cout << cost_adjlist[i][j].second << " ";
        }
        cout << endl;
    }
    cout << endl;

}

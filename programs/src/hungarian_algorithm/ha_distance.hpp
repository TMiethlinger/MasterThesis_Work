// Author:     Thomas Miethlinger B.Sc.
// Date:       21.10.2019
// Email:      t.miethlinger@protonmail.com
// University: Johannes Kepler University Linz, Austria

#include <cmath> // fabs
#include <vector> // std::vector

#pragma once

namespace ha_distance
{
    using std::vector;

    typedef vector<double> VD;

    // Distance functions

    double d_sum_3(VD &xi, VD &xj, double l_inv, double v_inv)
    {
        double dx, dy, dz, dvx, dvy, dvz;
        dx  = fabs(xi[0] - xj[0]);
        dy  = fabs(xi[1] - xj[1]);
        dz  = fabs(xi[2] - xj[2]);
        dvx = fabs(xi[3] - xj[3]);
        dvy = fabs(xi[4] - xj[4]);
        dvz = fabs(xi[5] - xj[5]);
        return l_inv * (dx*dx + dy*dy + dz*dz) + v_inv * (dvx*dvx + dvy*dvy + dvz*dvz);
    }

    double d_e_3(VD &xi, VD &xj, double l_inv, double v_inv)
    {
        double dx, dy, dz, drdl, dvx, dvy, dvz, dvdl;
        dx  = fabs(xi[0] - xj[0]);
        dy  = fabs(xi[1] - xj[1]);
        dz  = fabs(xi[2] - xj[2]);
        drdl = l_inv * (dx*dx + dy*dy + dz*dz);
        dvx = fabs(xi[3] - xj[3]);
        dvy = fabs(xi[4] - xj[4]);
        dvz = fabs(xi[5] - xj[5]);
        dvdl = v_inv * (dvx*dvx + dvy*dvy + dvz*dvz);
        return sqrt(drdl * drdl + dvdl * dvdl);
    }
};


    // TODO

/*  double euclidean_distance(std::vector<double> &ri, std::vector<double> &vi, std::vector<double> &rj, std::vector<double> &vj, double ratio);

    double euclidean_distance_position(std::vector<double> &ri, std::vector<double> &rj);
        
    double euclidean_distance_pbc(std::vector<double> &ri, std::vector<double> &vi, std::vector<double> &rj, std::vector<double> &vj, double L, double L_half);

    double euclidean_distance_position_pbc(std::vector<double> &ri, std::vector<double> &rj, double L, double L_half);*/

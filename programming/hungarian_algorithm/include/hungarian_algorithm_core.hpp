#include <algorithm>
#include <functional>
#include <utility>
#include <vector>
#include <list>
#include <set>
#include <cmath>
#include <numeric>
#include <limits>
#include <chrono>
#include <ctime>

#pragma once

class ParticleAssignment
{
    public:
        // Define public structs:
        struct SkeletonNode
        {
            SkeletonNode()
            {}
            SkeletonNode(double _cost, int _index_column) : cost(_cost), index_column(_index_column)
            {}
            double cost;
            int index_column;
        };
        struct CostObject
        {
            CostObject()
            {}
            CostObject(double _total_cost, std::vector<double>_cost_vector) : total_cost(_total_cost), cost_vector(_cost_vector)
            {}
            double total_cost;
            std::vector<double> cost_vector;
        };
        struct Result
        {
            Result(bool _success, ParticleAssignment::CostObject _cost, std::vector<int> &_q, std::vector<int> &_q_inv, double _time) :
            success(_success), cost(_cost), q(_q), q_inv(_q_inv), time(_time)
            {}
            bool success;
            ParticleAssignment::CostObject cost;
            std::vector<int> q;
            std::vector<int> q_inv;
            double time;
        };

        // Assignment methods
        static ParticleAssignment::Result MinCostMatching(std::vector<std::vector<double>> &Ri, std::vector<std::vector<double>> &Rj, int M, int B);
        static ParticleAssignment::Result MinCostMatchingMatrix(std::vector<std::vector<double>> &costmatrix);
        static ParticleAssignment::Result MinCostMatchingSkeleton(std::vector<std::vector<ParticleAssignment::SkeletonNode>> &costskeleton);

        // Creator methods
        // Solid wall BC
        static std::vector<std::vector<double>> CreateCostMatrix(std::vector<std::vector<double>> &Ri, std::vector<std::vector<double>> &Rj);
        static std::vector<std::vector<double>> CreateCostMatrix(std::vector<std::vector<double>> &Ri, std::vector<std::vector<double>> &Rj, int M);
        static std::vector<std::vector<ParticleAssignment::SkeletonNode>> CreateCostSkeleton(std::vector<std::vector<double>> &Ri, std::vector<std::vector<double>> &Rj, int M);

        // Periodic BC
        //static std::vector<std::vector<double>> CreateCostMatrixPBC(std::vector<std::vector<double>> &Ri, std::vector<std::vector<double>> &Rj, double L, double L2);

    private:

        // Compute cost
        static ParticleAssignment::CostObject CostFromMatchingMatrix(std::vector<std::vector<double>> &costmatrix, std::vector<int> &q);
        static ParticleAssignment::CostObject CostFromMatchingSkeleton(std::vector<std::vector<ParticleAssignment::SkeletonNode>> &costskeleton, std::vector<int> &q);
        static ParticleAssignment::CostObject CostFromMatchingVectors(std::vector<std::vector<double>> &Ri, std::vector<std::vector<double>> &Rj, std::vector<int> &q);

        // Metric
        static double EuclideanMetric(std::vector<double> &ri, std::vector<double> &rj);
		//static double EuclideanMetricPBC(std::vector<double> &ri, std::vector<double> &rj, double L, double L2);

        /* Object handling */
        // Private default constructor
        ParticleAssignment();
        // Private copy constructor
        ParticleAssignment(const ParticleAssignment& input);
        // Private copy assignment
        ParticleAssignment& operator=(const ParticleAssignment& input);

        // Precision which is used for all computations
        static constexpr double precision = 0.000000001;
        // Maximum cost, symbolizing infinity
        static constexpr double cost_max = (double)std::numeric_limits<int>::max();
};

#pragma once

#include "consts.h"
#include <Eigen/Dense>
#include <vector>
#include <list>

using std::vector;
using std::list;
using std::pair;
using Eigen::VectorXd;
using Eigen::MatrixXd;

void solve(double **Nodes, int P, int **Elements, int E, int **BoundElements, int Eb,
           vector<list<pair<int,int>>> ConjElements,
           VectorXd& Time, MatrixXd& T,
           MatrixXd& U, MatrixXd& sigma, MatrixXd& q);

void set_init_values(int P, MatrixXd& U, VectorXd& eps, MatrixXd& sigma, MatrixXd& omega,
                     MatrixXd& q, vector<PhaseType>& phase);

#ifndef VTKEXPORT_H
#define VTKEXPORT_H

#include <Eigen/Dense>
#include <vector>
#include <string>
#include "consts.h"

using std::vector;
using std::string;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Matrix3d;

void ExportTimeStep(std::string meshfile,  const string outfile, const VectorXd T,
                    const VectorXd& U, const VectorXd& sigma, const VectorXd& omega, const VectorXd& eps,
                     const VectorXd& q, vector<PhaseType> Phase, int P);


#endif // VTKEXPORT_H

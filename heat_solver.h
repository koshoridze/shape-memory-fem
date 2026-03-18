#ifndef HEAT_SOLVER_H
#define HEAT_SOLVER_H

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <list>
#include <vector>

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::Vector4d;
using Eigen::Vector3d;
using Eigen::Matrix2d;
using Eigen::Matrix3d;
using Eigen::Matrix4d;
using Eigen::MatrixXd;
using Eigen::PermutationMatrix;
using std::vector;
using std::list;

class SolverFEMHeat3d
{
public:
    SolverFEMHeat3d(double **Nodes, int P,
                      int **Elements, int E,
                      int **BoundElements, int Eb,
                      double k, double h,
                      const VectorXd& Time);
    void solve(MatrixXd& Temp);
private:
    void timestep(VectorXd& T, const VectorXd& Tpr, double tau);
    void applyBoundaryConditions(SpMat& K, VectorXd& b);
    void apply_Dirichlet(SpMat& A, VectorXd& b, int row, double value);
    void element_timematrix(int e, MatrixXd& Ce);
    void element_matrix(int e, MatrixXd& Ke);
    void element_vector(int e, MatrixXd& fe);
    void handle_element(int e);

    Matrix3d jacobian(double xi, double eta, double zeta);
    void calcCij(int i, int j, MatrixXd& Cij);
    void calcKij(int i, int j, MatrixXd& Kij);
    Matrix3d calcDl(double xi, double eta, double zeta);
    MatrixXd calcBil(int i, double xi, double eta, double zeta);

private:
    static constexpr int gn = 2;
    static constexpr double gauss[gn] = { -0.5773502691896257, 0.5773502691896257 };
    static constexpr double weights[gn] = {1, 1};

    /*static constexpr int gn = 3;
    static constexpr double gauss[gn] = { -0.7745966692414834, 0, 0.7745966692414834 };
    static constexpr double weights[gn] = {0.5555555555555556, 0.8888888888888888, 0.5555555555555556};*/

    static constexpr int n   = 8;
    static constexpr int dof = 1;
    double k, t;
    Matrix3d D;

    double **Nodes;
    int **Elements;
    int **BoundElements;
    int P, E;
    list<int> RemovedConstrains;

    vector<Vector3d> r0;

    vector<list<int>> ConjElements;

    VectorXd Time;
    double time_cur;
};

#endif // HEAT_SOLVER_H

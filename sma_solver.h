#pragma once

#include "consts.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <list>

typedef Eigen::SparseMatrix<double> SpMat;
using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::Matrix3d;
using Eigen::MatrixXd;
typedef Eigen::Triplet<double> Tr;
using std::vector;
using std::list;
using std::pair;

class SolverFEMSMA3d
{
public:
    SolverFEMSMA3d(double **Nodes, int P,
                    int **Elements, int E,
                    int **BoundElements, int Eb,
                    vector<list<pair<int,int>>> ConjElements,
                    const vector<PhaseType>& phasepr,
                    const VectorXd& Tpr, const VectorXd& T,
                    const VectorXd& sigmapr,
                    const VectorXd& epspr, const VectorXd& omegapr,
                    const VectorXd& qpr, int iter);

    void solve(VectorXd& Uk, VectorXd& Epsk, VectorXd& Sigmak,
               VectorXd& Omegak, VectorXd& Qk, vector<PhaseType>& Phasek);

private:
    void handle_node(int p);
    void handle_element(int e);
    void element_matrix(int e, MatrixXd& Ke);
    void element_vector(int e, MatrixXd& fe);

    void SMA_t(const VectorXd& sigma, const VectorXd& omega, double temp,
               double& tA, double& tM, double& T_sigma);
    void SMA_elastic(MatrixXd& CC, VectorXd& AC, double q);

    void SMA_direct(MatrixXd& CC, VectorXd& AC, const VectorXd& sigma,
                                     VectorXd& omega, VectorXd& domega,
                                     double temp, double q);

    void SMA_reverse(MatrixXd& CC, VectorXd& AC,
                                     const VectorXd& sigma, VectorXd& omega,
                                     double temp, double q);

    void calc_element_C(int e, double xi, double nu, double zeta, MatrixXd& CC, VectorXd& Ac);
    VectorXd calc_element_alpha(int e, double xi, double eta, double zeta);
    void calcKij(int e, int i, int j, MatrixXd& Kij);

    Matrix3d jacobian(double xi, double eta, double zeta);

    MatrixXd calcBi(int i, double xi, double eta, double zeta);

    void applyBoundaryConditions(SpMat& K, VectorXd& b);
    void apply_Dirichlet(SpMat& K, VectorXd& b, int row, double value);

    void construct_system(vector<Tr>& tripletList, VectorXd& b);

    void calc_strain(const VectorXd& Uk, VectorXd& Epsk, int p);
    void calc_stress(const VectorXd& Epsk, VectorXd& Sigmak, int p);
    void increment_parameters(const VectorXd& Sigmak, VectorXd& Omegak,
                               VectorXd& Qk, vector<PhaseType>& Phasek, int p);

    double stress_p_inv(const VectorXd& sigma);
    double stress_i_inv(const VectorXd& sigma);
    VectorXd stress_dev(const VectorXd& sigma, double sigma_p);

private:
    /*static constexpr int gn = 1;
    static constexpr double gauss[gn] = { 0 };
    static constexpr double weights[gn] = {2};*/

    /*static constexpr int gn = 2;
    double gauss[gn] = { -1/sqrt(3), 1/sqrt(3) };
    double weights[gn] = {1, 1};*/

    static constexpr int gn = 3;
    double gauss[gn] = { -sqrt(3./5.), 0, sqrt(3./5.) };
    double weights[gn] = {5./9., 8./9., 5./9.};

    /*static constexpr int gn = 4;
    double gauss[gn] = {-sqrt((3./7.)+(2./7.)*sqrt(6./5.)), -sqrt((3./7.)-(2./7.)*sqrt(6./5.)),
                         sqrt((3./7.)-(2./7.)*sqrt(6./5.)),  sqrt((3./7.)+(2./7.)*sqrt(6./5.))};
    double weights[gn] = {(18-sqrt(30))/36, (18+sqrt(30))/36,
                          (18+sqrt(30))/36, (18-sqrt(30))/36};*/

    static constexpr int n   = 8;
    static constexpr int dof = 3;
    static constexpr int ntens = 6;
    static constexpr int ndi = 3;

    double **Nodes;
    int **Elements;
    int **BoundElements;
    int P, E, Eb, iter;

    vector<Vector3d> r0;
    list<int> RemovedConstrains; /*!< Список с номерами узлов, к которым уже применили ГУ */
    vector<list<pair<int,int>>> ConjElements;

    VectorXd T, Tpr, Sigmapr, Epspr,
             Omegapr, Omega_med,
             Qpr;
    vector<PhaseType> Phasepr, Phasek;

    vector<MatrixXd> CC_n;
    vector<VectorXd> AC_n;
};

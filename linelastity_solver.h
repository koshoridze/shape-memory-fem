#pragma once

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <list>
#include <vector>

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> Trl;
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
using std::pair;

class SolverFEMElastity
{
public:
    SolverFEMElastity(double **Nodes, int P,
                      int **Elements, int E,
                      int **BoundElements, int Eb,
                      vector<list<pair<int,int>>> ConjElements,
                      double Em, double mu, double alpha);

    void solve(VectorXd& U, VectorXd& Eps, VectorXd& Sigma);

private:
    void applyBoundaryConditions(SpMat& K, VectorXd& b);
    void apply_Dirichlet(SpMat& A, VectorXd& b, int row, double value);

    void element_matrix(int e, MatrixXd& Ke);
    void element_vector(int e, MatrixXd& fe);
    void handle_element(int e);

    MatrixXd calcBi (int i, double xi, double eta, double zeta);
    Matrix3d jacobian(double xi, double eta, double zeta);
    void calcKij(int i, int j, MatrixXd& Kij);

    MatrixXd calcP(double x, double y, double z);
    void calc_strain(const VectorXd& Uk, VectorXd& Epsk, int p);
    void calc_stress(const VectorXd& Epsk, VectorXd& Sigmak, int p);

private:
    /*static constexpr int gn = 2; //число узлов интегрирования по Гауссу
    static constexpr double gauss[gn] = { -0.5773502691896257, 0.5773502691896257 };
    static constexpr double weights[gn] = {1, 1};*/

    static constexpr int gn = 3;
    static constexpr double gauss[gn] = { -0.7745966692414834, 0, 0.7745966692414834 };
    static constexpr double weights[gn] = {0.5555555555555556, 0.8888888888888888, 0.5555555555555556};


    static constexpr int n   = 8; /*!< Количество узлов в элементе */
    static constexpr int dof = 3; /*!< Количество степеней свободы в элементе */
    static constexpr int ntens = 6;
    double Em; /*!< Модуль Юнга */
    double nu; /*!< Коэффициент Пуассона */
    double G;
    MatrixXd D; /*!< Матрица жесткости для изотропного закона Гука, */

    double **Nodes; /*!< Массив всех узлов в сетке */
    int **Elements; /*!< Массив всех элементов в сетке */
    int **BoundElements; /*!< Массив граничных элементов */
    int P; /*!< Количество узлов */
    int E; /*!< Количество элементов */
    int Eb; /*!< Количество граничных элементов */
    vector<list<pair<int,int>>> ConjElements;
    list<int> RemovedConstrains; /*!< Список с номерами узлов, к которым уже применили ГУ */
    VectorXd Alpha;

    vector<Vector3d> r0; /*!< Массив из n радиус-векторов узлов в текущем элементе */
};


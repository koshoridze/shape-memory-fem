#include "heat_solver.h"
#include "shapefun.h"
#include "consts.h"
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>
#include <iostream>


/* \dot{T} = k $Delta T */
SolverFEMHeat3d::SolverFEMHeat3d(double **Nodes, int P,
                      int **Elements, int E,
                      int **BoundElements, int Eb,
                      double k, double h,
                      const VectorXd& Time)
{
    this->Nodes = Nodes;
    this->Elements = Elements;
    this->BoundElements = BoundElements;
    this->P  = P;
    this->E  = E;
    this->k = k;
    this->t  = h;

    r0.resize(n);

    D.setZero();
    D(0,0) = k;
    D(1,1) = k;
    D(2,2) = k;

    this->Time = Time;
}

void SolverFEMHeat3d::solve(MatrixXd& Temp)
{
    VectorXd Tpr(P), T(P);
    Tpr.setOnes();
    Tpr = T0 * Tpr;

    Temp.row(0) = Tpr;
    double tau1 = (T1-T0)/(N1-1);
    double tau2 = (N2 != 0) ? (T2-T1)/(N2) : 0.;
    double tau3 = (N3 != 0) ? (T3-T2)/(N3) : 0.;
    double tcur = T0;
    for (int i = 1; i < N1; i++)
    {
        tcur = (i != N1-1) ? tcur+tau1 : T1;
        Tpr.setOnes();
        Tpr = tcur * Tpr;
        Temp.row(i) = Tpr;
    }
    for (int i = N1; i < N1+N2; i++)
    {
        tcur = (i != N1+N2-1) ? tcur+tau2 : T2;
        Tpr.setOnes();
        Tpr = tcur * Tpr;
        Temp.row(i) = Tpr;
    }
    for (int i = N1+N2; i < NTime; i++)
    {
        tcur = (i != NTime-1) ? tcur+tau3 : T3;
        Tpr.setOnes();
        Tpr = tcur * Tpr;
        Temp.row(i) = Tpr;
    }


    for (int p = 0; p < P; p++)
    {
        /*if (abs(Nodes[p][0]) < 1e-5)
        {
            VectorXd d1(NTime);
            d1.setOnes();
            d1 = T0*d1;
            Temp.col(p) = d1;
        }*/
    }

    /*int Nt = Time.size();
    for (int i = 1; i < Nt; i++)
    {
        time_cur = Time(i);
        timestep(T, Tpr, Time(i)-Time(i-1));
        Temp.row(i) = T;
        Tpr = T;
    }*/
}

void SolverFEMHeat3d::timestep(VectorXd& Te, const VectorXd& Tpr, double tau)
{
    SpMat K(dof*P,dof*P);
    SpMat C(dof*P,dof*P);
    vector<T> tripletList, tripletList2;
    tripletList.reserve(3*P);
    tripletList2.reserve(2*P);
    VectorXd f(dof*P);
    f.setZero();
    MatrixXd Ke(dof*n, dof*n);
    MatrixXd Ce(dof*n, dof*n);
    MatrixXd fe(dof, n);

    for (int e = 0; e < E; e++)
    {
        handle_element(e);
        element_timematrix(e, Ce);
        element_matrix(e, Ke);
        element_vector(e, fe);
        for (int s = 0; s < n; s++)
        {
            for (int m = 0; m < n; m++)
            {
                for (int dx = 0; dx < dof; dx++)
                {
                    for (int dy = 0; dy < dof; dy++)
                    {
                        tripletList.push_back(
                                    T(dof*Elements[e][s]+dx, dof*Elements[e][m]+dy, Ke(dof*s+dx, dof*m+dy))
                                    );
                        tripletList2.push_back(
                                    T(dof*Elements[e][s]+dx, dof*Elements[e][m]+dy, Ce(dof*s+dx, dof*m+dy))
                                    );
                    }
                }
            }
            for (int dl = 0; dl < dof; dl++)
            {
                f(dof*Elements[e][s]+dl) += fe(dl,s);
            }
        }
    }
    K.setFromTriplets(tripletList.begin(), tripletList.end());
    C.setFromTriplets(tripletList2.begin(), tripletList2.end());

    /* Итоговая система */
    SpMat A = (1./tau)*C + K;
    VectorXd b = (1./tau)*C*Tpr + f;

    applyBoundaryConditions(A, b);

    A.makeCompressed();

    //Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int>> solver;
    //Eigen::BiCGSTAB<SpMat> solver;
    Eigen::ConjugateGradient<SpMat, Eigen::Lower|Eigen::Upper> solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Error: factorize(A)" << std::endl;
        return;
    }
    VectorXd Res = solver.solve(b);

    VectorXd r = A*Res - b;
    std::cout << "|r| =  " << r.norm() << std::endl;
    double maxabs = abs(r(0));
    for (int i = 1; i < P; i++)
    {
        if (maxabs < abs(r(i)))
            maxabs = abs(r(i));
    }
    std::cout << "max|r_i| = " << maxabs << std::endl;

    for (int i = 0; i < P; i++)
    {
        Te(i) = Res(dof*i);
    }

}

void SolverFEMHeat3d::applyBoundaryConditions(SpMat& A, VectorXd& b)
{
    const double delta = 1e-5;
    //const double l = 50;
    //double T0 = (time_cur < 1.5) ? 100 : -100;
    for (int p = 0; p < P; p++)
    {
        if (abs(Nodes[p][0]) < delta)
        {
            apply_Dirichlet(A, b, p, T0);
        }
        /*else if (abs(Nodes[p][0]-1) < delta)
        {
            apply_Dirichlet(A, b, p, 0);
        }*/
        /*else if (abs(Nodes[p][1]+l) < delta)
        {
            apply_Dirichlet(K, b, p, (Nodes[p][0]+l)/(2*l)+1);
        }
        else if (abs(Nodes[p][1]-l) < delta)
        {
            apply_Dirichlet(K, b, p, (Nodes[p][0]+l)/(2*l)+1);
        }*/
    }
    RemovedConstrains.clear();
}

void SolverFEMHeat3d::apply_Dirichlet(SpMat& K, VectorXd& b, int row, double value)
{
    if (std::find(RemovedConstrains.begin(), RemovedConstrains.end(), row) != RemovedConstrains.end())
        return;

    /* коррекция правой части */
    for (int i = 0; i < P; i++)
        if (i != row)
            b(i) -= K.coeff(i,row)*value;

    auto pr = [row](int i, int j, double)
    {
        return i != row && j != row;
    };
    K.prune(pr);
    K.coeffRef(row, row) = 1.;
    b(row) = value;
    RemovedConstrains.push_back(row);
}

void SolverFEMHeat3d::element_timematrix(int /*e*/, MatrixXd& Ce)
{
    MatrixXd Csm(dof, dof);
    Csm.setZero();

    for (int s = 0; s < n; s++)
    {
        for (int m = s; m < n; m++)
        {
            calcCij(s, m, Csm);
            for (int dx = 0; dx < dof; dx++)
            {
                for (int dy = 0; dy < dof; dy++)
                {
                    Ce(dof*s+dx, dof*m+dy) = Csm(dx,dy);
                    Ce(dof*m+dy, dof*s+dx) = Csm(dx,dy);
                }
            }
        }
    }
}

void SolverFEMHeat3d::element_matrix(int /*e*/, MatrixXd& Ke)
{
    MatrixXd Ksm(dof, dof);
    Ke.setZero();

    for (int s = 0; s < n; s++)
    {
        for (int m = s; m < n; m++)
        {
            calcKij(s, m, Ksm);
            for (int dx = 0; dx < dof; dx++)
            {
                for (int dy = 0; dy < dof; dy++)
                {
                    Ke(dof*s+dx, dof*m+dy) = Ksm(dx,dy);
                    Ke(dof*m+dy, dof*s+dx) = Ksm(dx,dy);
                }
            }
        }
    }
}

void SolverFEMHeat3d::element_vector(int e, MatrixXd& fe)
{
    fe.setZero();
    for (int i = 0; i < n; i++)
    {
        for (int p = 0; p < gn; p++)
            for (int q = 0; q < gn; q++)
                for (int l = 0; l < gn; l++)
                {
                    Vector3d r; r.setZero();
                    for (int j = 0; j < n; j++)
                    {
                        r += N(j, gauss[p], gauss[q], gauss[l])*r0[j];
                    }
                    fe(0,i) += weights[p]*weights[q]*weights[l]*
                            N(i,gauss[p],gauss[q], gauss[l]) *
                            abs(jacobian(gauss[p],gauss[q],gauss[l]).determinant())*
                            (-5);
                }
    }
}

void SolverFEMHeat3d::handle_element(int e)
{
    for (int p = 0; p < n; p++)
    {
        r0[p][0] = Nodes[Elements[e][p]][0];
        r0[p][1] = Nodes[Elements[e][p]][1];
        r0[p][2] = Nodes[Elements[e][p]][2];
    }

}

Matrix3d SolverFEMHeat3d::jacobian(double xi, double eta, double zeta)
{
    Matrix3d Je;
    Je.setZero();
    for (int i = 0; i < n; i++)
    {
        Je.row(0)  += dxiN  (i, xi, eta, zeta)*r0[i];
        Je.row(1)  += detaN (i, xi, eta, zeta)*r0[i];
        Je.row(2)  += dzetaN(i, xi, eta, zeta)*r0[i];
    }
    return Je;
}

void SolverFEMHeat3d::calcCij(int i, int j, MatrixXd& Cij)
{
    Cij.setZero();
    for (int p = 0; p < gn; p++)
       for (int q = 0; q < gn; q++)
          for (int r = 0; r < gn; r++)
          {
            Cij(0,0) += weights[p]*weights[q]*weights[r]*
                    abs(jacobian(gauss[p], gauss[q], gauss[r]).determinant())*
                    N(i, gauss[p], gauss[q], gauss[r]) * N(j, gauss[p], gauss[q], gauss[r]);
          }

}

void SolverFEMHeat3d::calcKij(int i, int j, MatrixXd& Kij)
{
    Kij.setZero();
    for (int p = 0; p < gn; p++)
       for (int q = 0; q < gn; q++)
          for (int r = 0; r < gn; r++)
          {
            Kij += weights[p]*weights[q]*weights[r]*
                    abs(jacobian(gauss[p], gauss[q], gauss[r]).determinant())*
                    calcBil(i,gauss[p], gauss[q], gauss[r]).transpose() *
                    calcDl(gauss[p], gauss[q], gauss[r]) *
                    calcBil(j,gauss[p], gauss[q], gauss[r]);
          }

}

MatrixXd SolverFEMHeat3d::calcBil(int i, double xi, double eta, double zeta)
{
    MatrixXd Bl(3,1);
    Bl(0,0) = dxiN  (i, xi, eta, zeta);
    Bl(1,0) = detaN (i, xi, eta, zeta);
    Bl(2,0) = dzetaN(i, xi, eta, zeta);

    return Bl;
}

Matrix3d SolverFEMHeat3d::calcDl(double xi, double eta, double zeta)
{
    Matrix3d invJe = jacobian(xi, eta, zeta).inverse();
    return (invJe.transpose()) * D * invJe;
}

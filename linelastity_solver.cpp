#include "linelastity_solver.h"
#include "shapefun.h"
#include "consts.h"
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>
#include <iostream>

SolverFEMElastity::SolverFEMElastity(double **Nodes, int P,
                      int **Elements, int E,
                      int **BoundElements, int Eb,
                      vector<list<pair<int,int>>> ConjElements,
                      double Em, double nu, double alpha)
{
    this->Nodes = Nodes;
    this->Elements = Elements;
    this->BoundElements = BoundElements;
    this->P  = P;
    this->E  = E;
    this->Eb = Eb;
    this->E  = E;
    this->Em = Em;
    this->nu = nu;
    this->ConjElements = ConjElements;

    r0.resize(n);

    Alpha.resize(ntens);
    Alpha.setZero();
    Alpha(0) = alpha;
    Alpha(1) = alpha;
    Alpha(2) = alpha;

    D.resize(6,6);
    D.setZero();
    double k1 = Em / ((1+nu)*(1-2*nu));
    D(0,0) = 1-nu;
    D(0,1) = nu; D(1,0) = D(0,1);
    D(0,2) = nu; D(2,0) = D(0,2);
    D(1,1) = 1-nu;
    D(1,2) = nu; D(2,1) = D(1,2);
    D(2,2) = 1-nu;
    D(3,3) = (1-2*nu)/2;
    D(4,4) = D(3,3);
    D(5,5) = D(3,3);
    D = k1 * D;
}

void SolverFEMElastity::solve(VectorXd& U, VectorXd& Eps, VectorXd& Sigma)
{
    // Numeration: u, v, w, th1, th2
    SpMat K(dof*P,dof*P);
    std::vector<Trl> tripletList;
    tripletList.reserve(4*P);
    VectorXd b(dof*P);
    b.setZero();
    MatrixXd Ke(dof*n, dof*n);
    MatrixXd fe(dof, n);

    for (int e = 0; e < E; e++)
    {
        handle_element(e);
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
                                    Trl(dof*Elements[e][s]+dx, dof*Elements[e][m]+dy, Ke(dof*s+dx, dof*m+dy))
                                    );
                    }
                }
            }
            for (int dl = 0; dl < dof; dl++)
            {
                b(dof*Elements[e][s]+dl) += fe(dl,s);
            }
        }
    }
    K.setFromTriplets(tripletList.begin(), tripletList.end());

    applyBoundaryConditions(K, b);
    for (int p = 0; p < P; p++)
    {
        if (abs(Nodes[p][2] - 2.) < 1e-4)
            b(dof*p+2) += Fz;
    }

    K.makeCompressed();

    //Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int>> solver;
    //Eigen::BiCGSTAB<SpMat> solver;
    Eigen::ConjugateGradient<SpMat, Eigen::Lower|Eigen::Upper> solver;
    solver.analyzePattern(K);
    solver.factorize(K);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Error: factorize(K)" << std::endl;
        return;
    }
    VectorXd Res = solver.solve(b);

    VectorXd r = K*Res - b;
    std::cout << "|r| =  " << r.norm() << std::endl;
    double maxabs = abs(r(0));
    for (int i = 1; i < P; i++)
    {
        if (maxabs < abs(r(i)))
            maxabs = abs(r(i));
    }
    std::cout << "max|r_i| = " << maxabs << std::endl;

    U = Res;

    for (int p = 0; p < P; p++)
        calc_strain(U, Eps, p);

    for (int p = 0; p < P; p++)
        calc_stress(Eps, Sigma, p);

}

void SolverFEMElastity::applyBoundaryConditions(SpMat& K, VectorXd& b)
{
    const double delta = 1e-6;
    //const double l = 25;
    for (int p = 0; p < P; p++)
    {
        if ( (abs(Nodes[p][0]-1) < delta) && (abs(Nodes[p][1]) < delta) && (abs(Nodes[p][2]-1.) < delta) )
        {
            apply_Dirichlet(K, b, dof*p  , 0);
            apply_Dirichlet(K, b, dof*p+1, 0);
            apply_Dirichlet(K, b, dof*p+2, 0);
        }
        else if ( (abs(Nodes[p][0]-1) < delta) && (abs(Nodes[p][1]-1) < delta) && (abs(Nodes[p][2]-1.) < delta) )
        {
            apply_Dirichlet(K, b, dof*p  , 0);
            apply_Dirichlet(K, b, dof*p+2, 0);
        }
        else if (abs(Nodes[p][2]-1.) < delta)
        {
            apply_Dirichlet(K, b, dof*p+2, 0);
        }
        if (abs(Nodes[p][2]-1.) < delta)
        {
            apply_Dirichlet(K, b, dof*p+0, 0);
            apply_Dirichlet(K, b, dof*p+1, 0);
            apply_Dirichlet(K, b, dof*p+2, 0);
        }
    }
}

void SolverFEMElastity::apply_Dirichlet(SpMat& K, VectorXd& b, int row, double value)
{
    if (std::find(RemovedConstrains.begin(), RemovedConstrains.end(), row) != RemovedConstrains.end())
        return;

    /* коррекция правой части */
    for (int i = 0; i < dof*P; i++)
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

void SolverFEMElastity::element_matrix(int /*e*/, MatrixXd& Ke)
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

void SolverFEMElastity::element_vector(int e, MatrixXd& fe)
{
    fe.setZero();

    for (int i = 0; i < n; i++)
    {
        //int p = Elements[e][i];
        for (int p = 0; p < gn; p++)
            for (int q = 0; q < gn; q++)
                for (int l = 0; l < gn; l++)
            {
                /*fe.col(i) += weights[p]*weights[q]*weights[l]*
                        abs(jacobian(gauss[p],gauss[q],gauss[l]).determinant())*
                        calcBi(i, gauss[p],gauss[q],gauss[l]).transpose() *
                        (dt*D*Alpha);*/
                fe(2,i) += weights[p]*weights[q]*weights[l]*
                                            N(i,gauss[p],gauss[q],gauss[l])*
                                            abs(jacobian(gauss[p],gauss[q],gauss[l]).determinant())*
                                            (-bz0);
            }
    }
}

void SolverFEMElastity::handle_element(int e)
{
    for (int p = 0; p < n; p++)
    {
        r0[p][0] = Nodes[Elements[e][p]][0];
        r0[p][1] = Nodes[Elements[e][p]][1];
        r0[p][2] = Nodes[Elements[e][p]][2];
    }
}

Matrix3d SolverFEMElastity::jacobian(double xi, double eta, double zeta)
{
    Matrix3d Je;
    Je.setZero();
    for (int i = 0; i < n; i++)
    {
        Je.row(0) += dxiN  (i, xi, eta, zeta)*r0[i];
        Je.row(1) += detaN (i, xi, eta, zeta)*r0[i];
        Je.row(2) += dzetaN(i, xi, eta, zeta)*r0[i];
    }
    return Je;
}

void SolverFEMElastity::calcKij(int i, int j, MatrixXd& Kij)
{
    Kij.setZero();
    for (int p = 0; p < gn; p++)
       for (int q = 0; q < gn; q++)
          for (int r = 0; r < gn; r++)
          {
            Kij += weights[p]*weights[q]*weights[r]*
                    abs(jacobian(gauss[p], gauss[q], gauss[r]).determinant())*
                    calcBi(i, gauss[p], gauss[q], gauss[r]).transpose() * D *
                    calcBi(j, gauss[p], gauss[q], gauss[r]);
          }

}

MatrixXd SolverFEMElastity::calcBi(int i, double xi, double eta, double zeta)
{
    MatrixXd B(6,3); B.setZero();
    Matrix3d invJe = jacobian(xi, eta, zeta).inverse();
    Vector3d dlN;
    dlN(0) = dxiN  (i, xi, eta, zeta);
    dlN(1) = detaN  (i, xi, eta, zeta);
    dlN(2) = dzetaN  (i, xi, eta, zeta);
    Vector3d dN = invJe * dlN;

    B(0,0) = dN(0);
    B(1,1) = dN(1);
    B(2,2) = dN(2);
    B(3,0) = B(1,1); B(3,1) = B(0,0);
    B(4,0) = B(2,2); B(4,2) = B(0,0);
    B(5,1) = B(2,2); B(5,2) = B(1,1);

    return B;
}

MatrixXd SolverFEMElastity::calcP(double x, double y, double z)
{
    MatrixXd P(1, 7);
    P(0,0) = 1;
    P(0,1) = x;
    P(0,2) = y;
    P(0,3) = z;
    P(0,4) = x*y;
    P(0,5) = x*z;
    P(0,6) = y*z;

    return P;
}

void SolverFEMElastity::calc_strain(const VectorXd& Uk, VectorXd& Epsk, int p)
{
    int ne = ConjElements[p].size();
    MatrixXd A(7,7); A.setZero();
    VectorXd b11(7), b22(7), b33(7),
        b12(7), b13(7), b23(7);
    b11.setZero(); b22.setZero();
    b12.setZero(); b13.setZero();
    b23.setZero(); b33.setZero();
    auto iter = ConjElements[p].begin();

    for (int ei = 0; ei < ne; ei++)
    {
        int e = (*iter).first; iter++;
        //int curi = (*iter).second; iter++;
        handle_element(e);
        for (int k = 0; k < gn; k++)
                        for (int q = 0; q < gn; q++)
                            for (int l = 0; l < gn; l++)
                            {
        double E11(0), E22(0), E33(0),
            E12(0), E13(0), E23(0);
        Vector3d r; r.setZero();

                    for (int j = 0; j < n; j++)
                    {
                        int p0 = Elements[e][j];
                        VectorXd ae(dof);
                        ae(0) = Uk(dof*p0  );
                        ae(1) = Uk(dof*p0+1);
                        ae(2) = Uk(dof*p0+2);;
                        VectorXd epsv = calcBi(j, gauss[k],gauss[q],gauss[l]) * ae;
                        E11 += epsv(0);
                        E22 += epsv(1);
                        E33 += epsv(2);
                        E12 += epsv(3);
                        E13 += epsv(4);
                        E23 += epsv(5);
                        r += N(j, gauss[k],gauss[q],gauss[l])*(r0[j]);
                    }
                    A += calcP(r(0), r(1), r(2)).transpose() * calcP(r(0), r(1), r(2));
                    b11 += E11*calcP(r(0), r(1), r(2)).transpose();
                    b22 += E22*calcP(r(0), r(1), r(2)).transpose();
                    b33 += E33*calcP(r(0), r(1), r(2)).transpose();
                    b12 += E12*calcP(r(0), r(1), r(2)).transpose();
                    b13 += E13*calcP(r(0), r(1), r(2)).transpose();
                    b23 += E23*calcP(r(0), r(1), r(2)).transpose();
                            }
    }
    MatrixXd invA = A.inverse();
    VectorXd P0 = calcP(Nodes[p][0], Nodes[p][1], Nodes[p][2]).transpose();
    Epsk(ntens*p+0) = P0.dot(invA * b11);
    Epsk(ntens*p+1) = P0.dot(invA * b22);
    Epsk(ntens*p+2) = P0.dot(invA * b33);
    Epsk(ntens*p+3) = P0.dot(invA * b12);
    Epsk(ntens*p+4) = P0.dot(invA * b13);
    Epsk(ntens*p+5) = P0.dot(invA * b23);

    /*int ne = ConjElements[p].size();
    auto iter = ConjElements[p].begin();
    VectorXd eps(ntens); eps.setZero();
    for (int ei = 0; ei < ne; ei++)
    {
        int e = (*iter).first;
        int curi = (*iter).second;
        iter++;

        for (int i = 0; i < n; i++)
        {
            int p0 = Elements[e][i];
            VectorXd ul(dof);
            ul(0) = Uk(dof*p0  );
            ul(1) = Uk(dof*p0+1);
            ul(2) = Uk(dof*p0+2);
            eps = eps + calcBi(i, xi_i(curi), eta_i(curi), eta_i(curi)) * ul;
        }
    }
    eps = eps / static_cast<double>(ne);

    for (int i = 0; i < ntens; i++)
        Epsk(ntens*p + i) = eps(i);*/
}

void SolverFEMElastity::calc_stress(const VectorXd& Epsk, VectorXd& Sigmak, int p)
{
    VectorXd strain(ntens), stress(ntens);
    for (int i = 0; i < ntens; i++)
    {
        strain(i) = Epsk(ntens*p+i);
    }
    stress = D*strain;
    for (int i = 0; i < ntens; i++)
        Sigmak(ntens*p+i) = stress(i);
}

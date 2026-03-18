#include "sma_solver.h"
#include "shapefun.h"
#include <vector>
#include <iostream>


SolverFEMSMA3d::SolverFEMSMA3d(double **Nodes, int P,
                int **Elements, int E,
                int **BoundElements, int Eb,
                vector<list<pair<int,int>>> ConjElements,
                const vector<PhaseType>& phasepr,
                const VectorXd& Tpr, const VectorXd& T,
                const VectorXd& sigmapr,
                const VectorXd& epspr, const VectorXd& omegapr,
                const VectorXd& Qpr, int iter)
{
    this->Nodes = Nodes;
    this->Elements = Elements;
    this->BoundElements = BoundElements;
    this->P  = P;
    this->E  = E;
    this->Eb = Eb;
    this->ConjElements = ConjElements;
    this->iter = iter;

    r0.resize(n);
    CC_n.resize(P);
    AC_n.resize(P);

    this->Tpr = Tpr;
    this->T   = T;

    this->Epspr = epspr;
    this->Sigmapr = sigmapr;
    this->Omegapr = omegapr;
    this->Qpr = Qpr;
    this->Phasepr = phasepr;

    Omega_med.resize(ntens*P);
}


void SolverFEMSMA3d::solve(VectorXd& Uk, VectorXd& Epsk, VectorXd& Sigmak,
                           VectorXd& Omegak, VectorXd& Qk, vector<PhaseType>& Phasek)
{
    SpMat K(dof*P, dof*P);
    VectorXd b(dof*P);
    b.setZero();
    vector<Tr> tripletList;
    tripletList.reserve(4*P);

    /* Calc СС, AC in each node */
    for (int p = 0; p < P; p++)
    {
        handle_node(p);
    }

    construct_system(tripletList, b);
    K.setFromTriplets(tripletList.begin(), tripletList.end());

    applyBoundaryConditions(K, b);
    for (int p = 0; p < P; p++)
    {
        if (abs(Nodes[p][2] - 2.) < 1e-4)
        {
            double F;
            /*if (iter <= N1)
                F = (static_cast<double>(iter)/N1)*Fz;
            else if (iter <= N1+N2)
            {
                F = (Fz/N2)*static_cast<double>(-iter+N1+N2);
            }
            else
                F = 0;*/
            F = Fz;
            b(dof*p+2) += F;
            //std::cout << F << std::endl;
        }
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
    Uk = solver.solve(b);
    VectorXd r = K*Uk - b;
    std::cout << "|r| =  " << r.norm() << std::endl;
    double maxabs = abs(r(0));
    for (int i = 1; i < P; i++)
    {
        if (maxabs < abs(r(i)))
            maxabs = abs(r(i));
    }
    std::cout << "max|r_i| = " << maxabs << std::endl;

    /* Calc new strains from new displacements */
    for (int p = 0; p < P; p++)
    {
        calc_strain(Uk, Epsk, p);
        calc_stress(Epsk, Sigmak, p);
        increment_parameters(Sigmak, Omegak, Qk, Phasek, p);
    }
}

void SolverFEMSMA3d::increment_parameters(const VectorXd& Sigmak, VectorXd& Omegak,
                                          VectorXd& Qk, vector<PhaseType>& Phasek, int p)
{
    PhaseType phase = Phasepr[p];
    double q = Qpr(p);
    double dq = 0;
    VectorXd sigmapr(ntens), sigma(ntens), omegapr(ntens), omega(ntens), domega(ntens);
    domega.setZero();
    for (int i = 0; i < ntens; i++)
    {
        sigmapr(i)   = Sigmapr(ntens*p+i);
        sigma(i)     = Sigmak(ntens*p+i);
        omegapr(i)   = Omegapr(ntens*p+i);
        omega(i)     = Omegapr(ntens*p+i); //Omega_med(ntens*p+i);
    }
    double temp = Tpr(p);
    double dtemp = T(p) - Tpr(p);

    double sigma_p       = stress_p_inv(sigma);
    double sigmapr_p     = stress_p_inv(sigmapr);
    double sigma_i       = stress_i_inv(sigma);
    double sigmapr_i     = stress_i_inv(sigmapr);
    VectorXd sigma_dev   = stress_dev(sigma, sigma_p);
    VectorXd sigmapr_dev = stress_dev(sigmapr, sigmapr_p);

    double T_sigma, tA, tA_pr, tM, tM_pr;
    SMA_t(sigmapr, omegapr, temp, tA_pr, tM_pr, T_sigma);

    // Direct martensitic transformation
    if (phase == DirectPhase)
    {
        if ((tM_pr >= 0) && (tM_pr <= 1))
        {
        double dF = (2/(sigma_0*sqrt(M_PI))) * exp(-pow(sigma_i/sigma_0,2));
        //double dF = erf(sigma_i/sigma_0);
        double tt = (sigma_i != 0.) ? 1.5*rhoD*dF / sigma_i : 0.;
        for (int i = 0; i < ntens; i++)
            domega(i) = tt*sigma_dev(i)*(sigma_i - sigmapr_i); // (2.14)
        }
    }

    omega += domega;
    SMA_t(sigma, omega, temp+dtemp, tA, tM, T_sigma);

    if (phase == DirectPhase)
    {
        if (tM <= 0)
            dq = -q;
        else if (tM >= 1) {
            dq = 1-q;
        }
        else
            dq = 0.5*(1-cos(M_PI*tM))-q;
    }

    if (phase == ReversePhase)
    {
        if (tA <= 0)
            dq = -q;
        else if (tA >= 1)
            dq = 1-q;
        else
            dq = 0.5*(1-cos(M_PI*tA))-q;
    }

    if ( (phase == DirectPhase) &&
        ( (tM >= 1)  && (tM_pr >= 1.) && (tM < tM_pr) ) )
        phase = ReversePhase;

    if ( (phase == ReversePhase) &&
        ( (tA <= 0) && (tA_pr <= 0.) && (tA > tA_pr) ) )
        phase = DirectPhase;

    Phasek[p] = phase;
    Qk(p) = q+dq;
    for (int i = 0; i < ntens; i++)
    {
        Omegak(ntens*p+i) = omega(i);
    }
}

void SolverFEMSMA3d::handle_node(int p)
{
    PhaseType phase = Phasepr[p];
    double q = Qpr(p);
    VectorXd sigmapr(ntens), omegapr(ntens), epspr(ntens);
    for (int i = 0; i < ntens; i++)
    {
        epspr(i)   = Epspr(ntens*p + i);
        sigmapr(i) = Sigmapr(ntens*p + i);
        omegapr(i) = Omegapr(ntens*p + i);
    }

    double sigma_p = stress_p_inv(sigmapr);
    double sigma_i = stress_i_inv(sigmapr);
    VectorXd sigma_dev = stress_dev(sigmapr, sigma_p);

    double beta = (sigma_i != 0) ? 1.5*rhoD*erf(sigma_i/sigma_0)/sigma_i : 0.;
    if (q == 0.)
        omegapr = beta*sigma_dev; // (2.9)

    double temp  = Tpr(p);
    double tA_pr, tM_pr, T_sigma;
    SMA_t(sigmapr, omegapr, temp, tA_pr, tM_pr, T_sigma);

    MatrixXd CC(ntens,ntens);
    VectorXd AC(ntens), domega(ntens);
    domega.setZero();

    //Direct martensitic transformation
    if (phase == DirectPhase)
    {
        if ((tM_pr < 0) || (tM_pr > 1))
        {
            SMA_elastic(CC, AC, q);
        }
        else
        {
            SMA_direct(CC, AC, sigmapr, omegapr, domega, temp,  q);
            q = 0.5*(1-cos(M_PI*tM_pr));
        }
    }

    // Reverse martensitic transformation
    if (phase == ReversePhase)
    {
        if ((tA_pr < 0) || (tA_pr > 1))
        {
            SMA_elastic(CC, AC, q);
        }
        else
        {
            SMA_reverse(CC, AC, sigmapr, omegapr, temp,  q);
            q = 0.5*(1-cos(M_PI*tA_pr));
        }
    }

    CC_n[p] = CC;
    AC_n[p] = AC;
    Qpr[p] = q;
}

void SolverFEMSMA3d::SMA_t(const VectorXd& sigma, const VectorXd& omega, double temp,
           double& tA, double& tM, double& T_sigma)
{
    double sigma_p = stress_p_inv(sigma);
    double sigma_i = stress_i_inv(sigma);
    VectorXd sigma_dev = stress_dev(sigma, sigma_p);

    double Z_sigma = (pow(sigma_p,2)*(1./Fm - 1./Fa) +
                      pow(sigma_i,2)*(1./Gm - 1./Ga)) / 6.;

    for (int i = 0; i < ntens; i++)
        Z_sigma += omega(i)*sigma_dev(i); //(2.5)

    T_sigma = temp - (Z_sigma + sigma_p*eps0) / dS; // (2.4)
    tA = 1. - (As - T_sigma) / (As - Af); // reverse transformation
    tM = (Ms - T_sigma) / (Ms - Mf); // direct transformation
}

void SolverFEMSMA3d::SMA_elastic(MatrixXd& CC, VectorXd& AC, double q)
{
    double Fk = 1/((q/Fm) + ((1-q)/Fa));
    double Gk = 1/((q/Gm) + ((1-q)/Ga));
    double alpha = q*alphaM + (1-q)*alphaA;

    AC.setZero();
    CC.setZero();

    for (int i = 0; i < ndi; i++) // (2.23)
    {
        AC(i) = alpha;
        for (int j = 0; j < ndi; j++)
            CC(i,j) = (2*Gk-Fk)/(6*Gk*Fk);
        CC(i,i) = (Fk+Gk)/(3*Fk*Gk);
    }
    for (int i = ndi; i < ntens; i++)
        CC(i,i) = 1./(2*Gk);

    CC = CC.inverse();

    /*for (int i = ndi; i < ntens; i++)
    {
        CC.col(i) *= 0.5;
        AC(i) *= 2;
    }*/
}

void SolverFEMSMA3d::SMA_direct(MatrixXd& CC, VectorXd& AC, const VectorXd& sigma,
                                 VectorXd& omega, VectorXd& domega,
                                 double temp, double q)
{
    double Fk = 1/((q/Fm) + ((1-q)/Fa));
    double Gk = 1/((q/Gm) + ((1-q)/Ga));
    double alpha = q*alphaM + (1-q)*alphaA;

    AC.setZero();
    CC.setZero();

    for (int i = 0; i < ndi; i++) // (2.23)
    {
        AC(i) = alpha;
        for (int j = 0; j < ndi; j++)
            CC(i,j) = (2*Gk-Fk)/(6*Gk*Fk);
        CC(i,i) = (Fk+Gk)/(3*Fk*Gk);
    }
    for (int i = ndi; i < ntens; i++)
        CC(i,i) = 1./(2*Gk);

    double sigma_p = stress_p_inv(sigma);
    double sigma_i = stress_i_inv(sigma);
    VectorXd sigma_dev = stress_dev(sigma, sigma_p);

    domega = omega;

    double GG = 0.5*((1/Gm) - (1/Ga));
    double eps0_pp = eps0 + (sigma_p/3)*((1/Fm) - (1/Fa));

    for (int i = 0; i < ntens; i++) // (2.24)
        omega(i) += GG*sigma_dev(i);
    for (int i = 0; i < ndi; i++)
        omega(i) += eps0_pp;

    //double dF = 2*exp(-pow(sigma_i,2)) / sqrt(M_PI);
    //dF = 0; // ??????

    double dF = (2/(sigma_0*sqrt(M_PI))) * exp(-pow(sigma_i/sigma_0,2));
    double tt = (alphaM - alphaA)*temp;
    VectorXd omega_s = omega;
    for (int i = 0; i < ndi; i++)
        omega_s(i) += tt; // (2.22)

    double fT = sqrt(abs(q*(1-q)))*M_PI / (Ms-Mf);
    double fN = fT / dS;
    double beta = 1.5*rhoD*dF;
    double delta = (sigma_i != 0.) ? 1.5*q*beta/pow(sigma_i, 2) : 0.;

    for (int i = 0; i < ntens; i++)
    {
        for (int j = 0; j < ntens; j++)
        {
            CC(i,j) += delta*(sigma_dev(i)*sigma_dev(j)) +
                        fN*omega_s(i)*(omega(j) + beta*sigma_dev(j)); // (2.21)
        }
    }

    for (int i = 0; i < ntens; i++)
        AC(i) -= fT*omega_s(i); // (2.22)

    CC = CC.inverse();

    /*for (int i = ndi; i < ntens; i++)
    {
        CC.col(i) *= 0.5;
        AC(i) *= 2;
    }*/

    omega = domega;
}

void SolverFEMSMA3d::SMA_reverse(MatrixXd& CC, VectorXd& AC,
                                 const VectorXd& sigma, VectorXd& omega,
                                 double temp, double q)
{
    double Fk = 1/((q/Fm) + ((1-q)/Fa));
    double Gk = 1/((q/Gm) + ((1-q)/Ga));
    double alpha = q*alphaM + (1-q)*alphaA;

    AC.setZero();
    CC.setZero();

    for (int i = 0; i < ndi; i++) // (2.23)
    {
        AC(i) = alpha;
        for (int j = 0; j < ndi; j++)
            CC(i,j) = (2*Gk-Fk)/(6*Gk*Fk);
        CC(i,i) = (Fk+Gk)/(3*Fk*Gk);
    }
    for (int i = ndi; i < ntens; i++)
        CC(i,i) = 1./(2*Gk);

    double sigma_p = stress_p_inv(sigma);
    //double sigma_i = stress_i_inv(sigma);
    VectorXd sigma_dev = stress_dev(sigma, sigma_p);

    VectorXd omega_temp = omega;

    double GG = 0.5*((1/Gm) - (1/Ga));
    double eps0_pp = eps0 + (sigma_p/3)*((1/Fm) - (1/Fa));

    for (int i = 0; i < ntens; i++) // (2.24)
        omega(i) += GG*sigma_dev(i);
    for (int i = 0; i < ndi; i++)
        omega(i) += eps0_pp;

    double tt = (alphaM - alphaA)*temp;
    VectorXd omega_s = omega;
    for (int i = 0; i < ndi; i++)
        omega_s(i) += tt;

    double fT = sqrt(abs(q*(1-q)))*M_PI / (As-Af);
    double fN = fT / dS;

    for (int i = 0; i < ntens; i++)
    {
        for (int j = 0; j < ntens; j++)
        {
            CC(i,j) -= fN*omega_s(i)*omega(j);
        }
    }

    for (int i = 0; i < ntens; i++)
        AC(i) += fT*omega_s(i);

    //for (int i = ndi; i < ntens; i++)
    //    CC.row(i) = 2 * CC.row(i);

    CC = CC.inverse();

    omega = omega_temp;

    /*dq = -dtemp*dS;
    for (int i = 0; i < ntens; i++)
        dq += omega(i)*dsigma(i);
    dq = dq*fT / dS;*/
}


void SolverFEMSMA3d::construct_system(vector<Tr>& tripletList, VectorXd& b)
{
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
                                        Tr(dof*Elements[e][s]+dx, dof*Elements[e][m]+dy, Ke(dof*s+dx, dof*m+dy))
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
}

void SolverFEMSMA3d::element_matrix(int e, MatrixXd& Ke)
{
    MatrixXd Ksm(dof, dof);
    Ke.setZero();

    for (int s = 0; s < n; s++)
    {
        for (int m = s; m < n; m++)
        {
            calcKij(e, s, m, Ksm);
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

void SolverFEMSMA3d::element_vector(int e, MatrixXd& fe)
{
    fe.setZero();

    MatrixXd C(ntens, ntens);
    VectorXd stress0(ntens), strain0(ntens), AC(ntens);
    double dt;
    for (int i = 0; i < n; i++)
    {
        for (int p = 0; p < gn; p++)
            for (int q = 0; q < gn; q++)
                for (int l = 0; l < gn; l++)
            {
                dt = T(Elements[e][i]) - Tpr(Elements[e][i]);
                calc_element_C(e, gauss[p], gauss[q], gauss[l], C, AC);
                //C = calc_element_C(e, gauss[p], gauss[q], gauss[l]);
                //AC = calc_element_alpha(e, gauss[p], gauss[q], gauss[l]);
                for (int k = 0; k < ntens; k++)
                {
                    strain0(k) = Epspr  (ntens*Elements[e][i]+k);
                    stress0(k) = Sigmapr(ntens*Elements[e][i]+k);
                }
                fe.col(i) += weights[p]*weights[q]*weights[l]*
                        abs(jacobian(gauss[p],gauss[q],gauss[l]).determinant())*
                        calcBi(i, gauss[p],gauss[q],gauss[l]).transpose() *
                        (
                            -stress0 + C*strain0 + dt*C*AC
                        );
                fe(2,i) += weights[p]*weights[q]*weights[l]*
                                            N(i,gauss[p],gauss[q],gauss[l])*
                                            abs(jacobian(gauss[p],gauss[q],gauss[l]).determinant())*
                                            (-bz0);
            }
    }
}

void SolverFEMSMA3d::calc_element_C(int e, double xi, double eta, double zeta,
                                    MatrixXd& CC, VectorXd& AC)
{
    CC.setZero();
    AC.setZero();
    for (int i = 0; i < n; i++)
    {
        CC += N(i, xi, eta, zeta)*CC_n[Elements[e][i]];
        AC += N(i, xi, eta, zeta)*AC_n[Elements[e][i]];
    }

    /*PhaseType phase = Phasepr[Elements[e][0]];
    for (int k = 0; k < n; k++)
    {
        if (Phasepr[Elements[e][k]] != phase)
        {
            for (int i = 0; i < n; i++)
            {
                CC += N(i, xi, eta, zeta)*CC_n[Elements[e][i]];
                AC += N(i, xi, eta, zeta)*AC_n[Elements[e][i]];
            }
            return;
        }
    }

    double q = 0;
    double temp = 0;
    VectorXd sigmapr(ntens), omegapr(ntens), epspr(ntens);
    sigmapr.setZero(); omegapr.setZero(); epspr.setZero();
    for (int k = 0; k < n; k++)
    {
        int p = Elements[e][k];
        q += N(k, xi, eta, zeta)*Qpr(p);
        temp += N(k, xi, eta, zeta)*Tpr(p);
        for (int i = 0; i < ntens; i++)
        {
            epspr(i)   += N(k, xi, eta, zeta)*Epspr(ntens*p + i);
            sigmapr(i) += N(k, xi, eta, zeta)*Sigmapr(ntens*p + i);
            omegapr(i) += N(k, xi, eta, zeta)*Omegapr(ntens*p + i);
        }
    }

    double sigma_p = stress_p_inv(sigmapr);
    double sigma_i = stress_i_inv(sigmapr);
    VectorXd sigma_dev = stress_dev(sigmapr, sigma_p);

    double beta = (sigma_i != 0) ? 1.5*rhoD*erf(sigma_i/sigma_0)/sigma_i : 0.;
    if (q == 0.)
        omegapr = beta*sigma_dev; // (2.9)

    double tA_pr, tM_pr, T_sigma;
    SMA_t(sigmapr, omegapr, temp, tA_pr, tM_pr, T_sigma);

    VectorXd  domega(ntens);
    domega.setZero();
    // Прямое мартенситное превращение
    if (phase == DirectPhase)
    {
        if ((tM_pr < 0) || (tM_pr > 1))
        {
            SMA_elastic(CC, AC, q);
        }
        else
        {
            SMA_direct(CC, AC, sigmapr, omegapr, domega, temp,  q);
            //q = 0.5*(1-cos(M_PI*tM_pr));
        }
    }

    // Обратное мартенситное превращение
    if (phase == ReversePhase)
    {
        if ((tA_pr < 0) || (tA_pr > 1))
        {
            SMA_elastic(CC, AC, q);
        }
        else
        {
            SMA_reverse(CC, AC, sigmapr, omegapr, temp, q);
            //q = 0.5*(1-cos(M_PI*tA_pr));
        }
    }*/

}

VectorXd SolverFEMSMA3d::calc_element_alpha(int e, double xi, double eta, double zeta)
{
    VectorXd AC(ntens);
    AC.setZero();
    for (int i = 0; i < n; i++)
    {
        AC += N(i, xi, eta, zeta)*AC_n[Elements[e][i]];
    }
    return AC;
}

void SolverFEMSMA3d::calcKij(int e, int i, int j, MatrixXd& Kij)
{
    Kij.setZero();
    MatrixXd CC(ntens, ntens);
    VectorXd AC(ntens);
    for (int p = 0; p < gn; p++)
       for (int q = 0; q < gn; q++)
          for (int r = 0; r < gn; r++)
          {
            calc_element_C(e, gauss[p], gauss[q], gauss[r], CC, AC);
            Kij += weights[p]*weights[q]*weights[r]*
                    abs(jacobian(gauss[p], gauss[q], gauss[r]).determinant())*
                    calcBi(i, gauss[p], gauss[q], gauss[r]).transpose() *
                    CC *
                    calcBi(j, gauss[p], gauss[q], gauss[r]);
          }

}

void SolverFEMSMA3d::handle_element(int e)
{
    for (int p = 0; p < n; p++)
    {
        r0[p][0] = Nodes[Elements[e][p]][0];
        r0[p][1] = Nodes[Elements[e][p]][1];
        r0[p][2] = Nodes[Elements[e][p]][2];
    }
}

Matrix3d SolverFEMSMA3d::jacobian(double xi, double eta, double zeta)
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

MatrixXd SolverFEMSMA3d::calcBi(int i, double xi, double eta, double zeta)
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

void SolverFEMSMA3d::applyBoundaryConditions(SpMat& K, VectorXd& b)
{
    double delta = 1e-6;
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

        /*if (abs(Nodes[p][2]-1.) < delta)
        {
            apply_Dirichlet(K, b, dof*p+0, 0);
            apply_Dirichlet(K, b, dof*p+1, 0);
            apply_Dirichlet(K, b, dof*p+2, 0);
        }*/
    }
}

void SolverFEMSMA3d::apply_Dirichlet(SpMat& K, VectorXd& b, int row, double value)
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

MatrixXd calcP(double x, double y, double z)
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

void SolverFEMSMA3d::calc_strain(const VectorXd& Uk, VectorXd& Epsk, int p)
{
    // Number of elements connected to node p
    int ne = ConjElements[p].size();

     // Least-squares system for local strain reconstruction
    MatrixXd A(7,7); A.setZero();

    // Right-hand sides for individual strain components
    VectorXd b11(7), b22(7), b33(7),
        b12(7), b13(7), b23(7);
    b11.setZero(); b22.setZero();
    b12.setZero(); b13.setZero();
    b23.setZero(); b33.setZero();


    auto iter = ConjElements[p].begin();

    // Loop over all elements adjacent to node p
    for (int ei = 0; ei < ne; ei++)
    {
        int e = (*iter).first; iter++;
        //int curi = (*iter).second; iter++;
        handle_element(e);
        // Loop over Gauss integration points in the element
        for (int k = 0; k < gn; k++)
            for (int q = 0; q < gn; q++)
                for (int l = 0; l < gn; l++)
                {
                    // Strain tensor components at the current Gauss point
                    double E11(0), E22(0), E33(0),
                        E12(0), E13(0), E23(0);
                    Vector3d r; r.setZero();

                    // Accumulate strain from element nodal displacements
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
                    // Assemble least-squares matrix and right-hand sides
                    A += calcP(r(0), r(1), r(2)).transpose() * calcP(r(0), r(1), r(2));
                    b11 += E11*calcP(r(0), r(1), r(2)).transpose();
                    b22 += E22*calcP(r(0), r(1), r(2)).transpose();
                    b33 += E33*calcP(r(0), r(1), r(2)).transpose();
                    b12 += E12*calcP(r(0), r(1), r(2)).transpose();
                    b13 += E13*calcP(r(0), r(1), r(2)).transpose();
                    b23 += E23*calcP(r(0), r(1), r(2)).transpose();
                            }
    }
    // Solve for polynomial coefficients of the reconstructed strain field
    MatrixXd invA = A.inverse();
    // Evaluate reconstructed strain at node p
    VectorXd P0 = calcP(Nodes[p][0], Nodes[p][1], Nodes[p][2]).transpose();
    Epsk(ntens*p+0) = P0.dot(invA * b11);
    Epsk(ntens*p+1) = P0.dot(invA * b22);
    Epsk(ntens*p+2) = P0.dot(invA * b33);
    Epsk(ntens*p+3) = P0.dot(invA * b12);
    Epsk(ntens*p+4) = P0.dot(invA * b13);
    Epsk(ntens*p+5) = P0.dot(invA * b23);


    /* Simpler alternative: strain averaging over adjacent elements
     * int ne = ConjElements[p].size();
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

void SolverFEMSMA3d::calc_stress(const VectorXd& Epsk, VectorXd& Sigmak, int p)
{
    VectorXd strain(ntens), strain0(ntens), stress(ntens), stress0(ntens);
    for (int i = 0; i < ntens; i++)
    {
        strain0(i) = Epspr  (ntens*p+i);
        strain (i) = Epsk   (ntens*p+i);
        stress0(i) = Sigmapr(ntens*p+i);
    }
    double dt = T(p)-Tpr(p);
    stress = stress0 + CC_n[p]*( strain - strain0 - dt*AC_n[p] ); // (2.29)
    for (int i = 0; i < ntens; i++)
        Sigmak(ntens*p+i) = stress(i);
}

double SolverFEMSMA3d::stress_p_inv(const VectorXd& sigma)
{
    // First stress invariant: trace of the stress tensor
    double sigma_p = 0;
    for (int i = 0; i < ndi; i++)
        sigma_p += sigma(i);
    return sigma_p;
}

double SolverFEMSMA3d::stress_i_inv(const VectorXd& sigma)
{
    // Second invariant of the deviatoric stress tensor
    // (equivalent von Mises stress)
    double sigma_i = 0;
    double sigma_p = stress_p_inv(sigma);
    for (int i = 0; i < ndi; i++)
        sigma_i += pow(sigma(i)-(sigma_p/3),2);
    for (int i = ndi; i < ntens; i++)
        sigma_i += 2*pow(sigma(i),2);
    return sqrt(3*sigma_i/2);
}

VectorXd SolverFEMSMA3d::stress_dev(const VectorXd& sigma, double sigma_p)
{
    // Deviatoric part of the stress tensor
    VectorXd sigma_dev = sigma;
    for (int i = 0; i < ndi; i++)
        sigma_dev(i) -= sigma_p / 3.;
    return sigma_dev;
}

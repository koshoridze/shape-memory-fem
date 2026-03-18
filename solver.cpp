#include "solver.h"
#include "consts.h"
#include "heat_solver.h"
#include "linelastity_solver.h"
#include "sma_solver.h"
#include "vtkexport.h"
#include <iostream>

// Initialize all primary field variables and internal state
void set_init_values(int P, MatrixXd& U, VectorXd& eps, MatrixXd& sigma, VectorXd& omega,
                     MatrixXd& q, vector<PhaseType>& phase)
{
    // Set initial phase for all material points (assume forward transformation)
    for (int i = 0; i < P; i++)
        phase[i] = DirectPhase;

    VectorXd temp(3*P);
    temp.setZero();

    // Initial displacement field (t = 0)
    U.row(0) = temp;

    // Initialize strain, transformation strain (omega), and stress
    temp.resize(6*P);
    temp.setZero();
    eps = temp;          // total strain at previous step
    omega = temp;        // transformation strain
    sigma.row(0) = temp; // initial stress state

    // Initialize martensite volume fraction
    temp.resize(P);
    temp.setOnes();
    temp = q0 * temp;
    q.row(0) = temp;
}

// Main FEM solver routine for thermo-mechanical SMA problem
void solve(double **Nodes, int P, int **Elements, int E, int **BoundElements, int Eb,
           vector<list<pair<int,int>>> ConjElements,
           VectorXd& Time, MatrixXd& T,
           MatrixXd& U, MatrixXd& sigma, MatrixXd& q)
{
    // Path to mesh (used for exporting results)
    string meshfile = "/home/georgiy/Documents/Учёба/Магистратура/Дипломная работа/ShapeMemory3d/sma3d/ShapeMemory3d/mesh.vtk";

    // Allocate storage for time history of fields
    Time.resize(NTime);
    T.resize(NTime, P);         // temperature field
    U.resize(NTime, 3*P);       // displacement field
    q.resize(NTime, P);         // martensite fraction
    sigma.resize(NTime, 6*P);   // stress tensor

    // Previous-step strain and transformation strain
    VectorXd epspr(6*P), omegapr(6*P);
    omegapr.setZero();
    vector<PhaseType> phasepr(P), phasek(P);

    // Time discretization
    double tau = Time_end / static_cast<double>(NTime-1);
    for (int i = 0; i < NTime; i++)
        Time(i) = (i != NTime-1) ? i*tau : Time_end;

    // Solve heat conduction problem (temperature field evolution).
    // NOTE: In the current version, the temperature field is prescribed
    // as a spatially uniform, piecewise-linear function of time.
    // The actual FEM heat solver (time-stepping with assembly) is implemented
    // but currently disabled in the code.
    SolverFEMHeat3d aHeatSolver(Nodes, P, Elements, E,
                                BoundElements, Eb,
                                1, 1e-2, Time);
    aHeatSolver.solve(T);

    std::cout << std::endl << std::endl<< std::endl;
    std::cout << "===== END OF TEMPERATURE SOLUTION  =====" << std::endl;
    std::cout << std::endl << std::endl<< std::endl;

    // Initialize mechanical and internal variables
    set_init_values(P, U, epspr, sigma, omegapr, q, phasepr);

     // Current step variables
    VectorXd Uk(3*P), Qk(P), Sigmak(6*P);

    /*SolverFEMElastity aLinSolver(Nodes, P, Elements, E,
                                 BoundElements, Eb,
                                 ConjElements, Ea, nu,
                                 alphaA);
    aLinSolver.solve(Uk, epspr, Sigmak);
    U.row(0) = Uk;
    sigma.row(0) = Sigmak;*/

    // Optional export of initial state
    //string filename = "/home/georgiy/Documents/Учёба/Магистратура/Дипломная работа/ShapeMemory3d/sma3d/ShapeMemory3d/results/data"
    //        + std::to_string(0) + ".vtk";
    //ExportTimeStep(meshfile, filename, T.row(0), U.row(0), sigma.row(0), omegapr, epspr, q.row(0), phasepr, P);

    // Main loop
    for (int i = 1; i < NTime; i++)
    {
        // Construct SMA FEM solver for current time step
        SolverFEMSMA3d aSMASolver(Nodes, P, Elements, E, BoundElements, Eb,
                                  ConjElements, phasepr, T.row(i-1), T.row(i),
                                  sigma.row(i-1), epspr, omegapr, q.row(i-1), i-1);

        aSMASolver.solve(Uk, epspr, Sigmak, omegapr, Qk, phasepr);

        //std::cout << T(i,0) << std::endl;
        std::cout << i-1 << std::endl;
        std::cout << "phase = " << phasepr[0] << std::endl << std::endl;

        U.row(i) = Uk;
        q.row(i) = Qk;
        sigma.row(i) = Sigmak;

         // Export results for visualization (VTK)
        string filename = "/home/georgiy/Documents/Учёба/Магистратура/Дипломная работа/ShapeMemory3d/sma3d/ShapeMemory3d/results/data"
                + std::to_string(i-1) + ".vtk";
        ExportTimeStep(meshfile, filename, T.row(i), U.row(i), sigma.row(i), omegapr, epspr,  q.row(i), phasepr, P);
    }
}

#include "importmesh.h"
#include "solver.h"
#include "vtkexport.h"
#include <iostream>
#include <list>
#include <cmath>
#include <vector>
#include <Eigen/Dense>

void FindAllConjElements(int P, double **Nodes, int E, int **Elements,
                         std::vector<std::list<std::pair<int,int>>>& ConjElements)
{
    const double eps = 1e-5;
    for (int p = 0; p < P; p++)
    {
        for (int e = 0; e < E; e++)
        {
            for (int i = 0; i < 8; i++)
            {
                int p0 = Elements[e][i];
                if ( (abs(Nodes[p0][0] - Nodes[p][0]) < eps) &&
                     (abs(Nodes[p0][1] - Nodes[p][1]) < eps) &&
                     (abs(Nodes[p0][2] - Nodes[p][2]) < eps) )
                {
                    ConjElements[p].push_back(std::pair<int,int>(e, i));
                    break;
                }
            }
        }
    }
}

int main()
{
    /* Mesh data */
    int P, E, Eb;
    double **Nodes;
    int **BoundElements;
    int **Elements;

    ImportMesh("/home/georgiy/Documents/Учёба/Магистратура/Дипломная работа/ShapeMemory3d/sma3d/ShapeMemory3d/mesh.vtk", &Nodes, P, &Elements, E, &BoundElements, Eb);

    std::vector<std::list<std::pair<int,int>>> ConjElements(P);
    FindAllConjElements(P, Nodes, E, Elements, ConjElements);

    VectorXd Time;
    MatrixXd T, U, sigma, q;

    solve(Nodes, P, Elements, E, BoundElements, Eb, ConjElements,
          Time, T, U, sigma, q);

    //ExportVtk("/home/georgiy/Desktop/ShapeMemory3d/sma3d/ShapeMemory3d/mesh.vtk", Time, T, U, sigma, q, NTime, P);
}

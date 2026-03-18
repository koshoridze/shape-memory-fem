#include "importmesh.h"

#include <string>
#include <iostream>
#include <fstream>

int ImportMesh(std::string fileName, double ***Nodes, int &P,
                int ***Elements, int &E, int ***BoundElements, int &Ed)
{
    int i;
    std::fstream meshFile(fileName, std::ios_base::in);

    if (!meshFile.is_open())
    {
        std::cerr << "Error while opening mesh file" << std::endl;
        return 0;
    }

    std::string line, skeep_me;
    for (i = 0; i < 4; i++)
        getline(meshFile, line);
    meshFile >> skeep_me;
    meshFile >> P;
    getline(meshFile, line);

    *Nodes = (double**) malloc(P*sizeof(double*));
    for (i = 0; i < P; i++)
    {
        (*Nodes)[i] = (double *) malloc(3*sizeof(double));
        meshFile >> (*Nodes)[i][0] >> (*Nodes)[i][1] >> (*Nodes)[i][2];
        getline(meshFile, line);
    }
    getline(meshFile, line);
    meshFile >> skeep_me;
    meshFile >> E;
    getline(meshFile, line);
    auto pos = meshFile.tellg();
    int d; Ed = 0;
    while (1) // !!
    {
        meshFile >> d;
        if ((d == 1) || (d == 2) || (d == 3))
            pos = meshFile.tellg();
        if (d == 4)
            Ed++;
        else if (d == 8)
            break;
        E--;
        getline(meshFile, line);
    }

    meshFile.seekg(pos);
    getline(meshFile, line);
    *BoundElements = (int **) malloc(Ed*sizeof(int*));
    for (i = 0; i < Ed; i++)
    {
        (*BoundElements)[i] = (int *) malloc(4*sizeof(int));
        meshFile >> d >> (*BoundElements)[i][0] >> (*BoundElements)[i][1] >> (*BoundElements)[i][2] >> (*BoundElements)[i][3];
        getline(meshFile, line);
    }

    *Elements = (int**) malloc(E*sizeof(int*));
    for (i = 0; i < E; i++)
    {
        (*Elements)[i] = (int *) malloc(8*sizeof(int));
        meshFile  >> d >> (*Elements)[i][0] >> (*Elements)[i][1] >> (*Elements)[i][2] >> (*Elements)[i][3] >>
                (*Elements)[i][4] >> (*Elements)[i][5] >> (*Elements)[i][6] >> (*Elements)[i][7];
        getline(meshFile, line);
    }
    meshFile.close();

    return 1;
}

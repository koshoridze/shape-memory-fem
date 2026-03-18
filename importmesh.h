#ifndef IMPORTMESH_H
#define IMPORTMESH_H

#include <string>
#include <list>

int ImportMesh(std::string fileName, double ***Nodes, int &P,
                int ***Elements, int &E, int ***BoundElements, int &Ed);
#endif

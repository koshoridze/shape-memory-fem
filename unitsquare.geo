

lc = 0.1;
a = 10;
b = 1;
c = 1;


Point(1) = {0, 0, 0, lc};
Point(2) = {0, b,  0, lc};
Point(3) = {a, b, 0, lc};
Point(4) = {a, 0, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {4,3,2,1};

Plane Surface(1) = {1};

Extrude {0,0,c} { Surface{1};};

Physical Surface("My surface") = {2};
Physical Volume("My volume")  = {1};

Mesh.SaveAll = 1;
Mesh.ElementOrder=1;
#Mesh.Algorithm = 6;
Mesh.RecombineAll = 1;

Mesh 3;
Save "mesh.vtk";

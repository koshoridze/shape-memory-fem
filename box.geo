SetFactory("OpenCASCADE");
Lbox = 1;
Hbox = 1;

Wref = 1;
size_el = 0.6;
n_pts_big = Lbox/size_el + 1.;
n_pts_small = Wref/size_el + 1.;

Rectangle(1) = {0, 0, 1, Lbox, Hbox};

Transfinite Curve {:} = n_pts_big Using Progression 1;
Transfinite Surface {1} Alternate;

Extrude {0, 0, Wref} {
   Surface{:}; Layers{n_pts_small}; Recombine;
} 

Mesh.SaveAll = 1;
Mesh.RecombineAll = 1;

Mesh 3;
Save "mesh.vtk";

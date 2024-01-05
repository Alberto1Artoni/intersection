//Inputs
w=0.5;//m
h=0.5;//m
r=1;//m
gridsize=0.1;
gridsize2=0.1;

Point(1) = {-w, -h, 0, gridsize};
Point(2) = {w, -h, 0, gridsize};
Point(3) = {w, h, 0, gridsize};
Point(4) = {-w, h, 0, gridsize};
Point(5) = {0, 0, 0, gridsize2};
Point(6) = {r, 0, 0, gridsize2};
Point(7) = {0, r, 0, gridsize2};
Point(8) = {-r, 0, 0, gridsize2};
Point(9) = {0, -r, 0, gridsize2};
Line(1) = {4, 3};
Line(2) = {3, 2};
Line(3) = {2, 1};
Line(4) = {1, 4};
Line Loop(1) = {4, 1, 2, 3};
Plane Surface(1) = {1};

surfaceVector[] = Extrude {0, 0, 0.1} {
Surface{1};
Layers{1};
Recombine;};

Physical Volume("internal")=surfaceVector[1];
Physical Surface("left")=surfaceVector[2];
Physical Surface("right")=surfaceVector[4];
Physical Surface("bottom")=surfaceVector[5];
Physical Surface("top")=surfaceVector[3];
Physical Surface("back")={1};
Physical Surface("front")={surfaceVector[0]};

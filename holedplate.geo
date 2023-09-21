// Gmsh project created on Fri Aug  4 10:27:14 2023
SetFactory("OpenCASCADE");
a=.5;
Point(1) = {-1, -1, 0, a};
//+
Point(2) = {1, -1, 0, a};
//+
Point(3) = {1, 1, 0, a};
//+
Point(4) = {-1, 1, 0, a};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
Line(4) = {4, 1};
//+
Circle(10) = {0,0, 0., 0.2, 0, 2*Pi};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Curve Loop(2) = {10};
//+
Plane Surface(1) = {1, 2};
//+
Physical Curve("Right") = {2};
Physical Curve("Left") = {4};
Physical Curve("Top") = {3};
Physical Curve("Bottom") = {1};
Physical Curve("Circle") = {10};
Physical Surface("Domain") = {1};


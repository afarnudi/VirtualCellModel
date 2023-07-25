// Gmsh project created on Sat Aug 17 15:33:32 2019
SetFactory("OpenCASCADE");
//+
Point(1) = {50, 50, 0, 1.0};
//+
Point(2) = {-50, 50, 0, 1.0};
//+
Point(3) = {-50, -50, 0, 1.0};
//+
Point(4) = {50, -50, 0, 1.0};
//+
Line(1) = {2, 3};
//+
Line(2) = {2, 1};
//+
Line(3) = {1, 4};
//+
Line(4) = {4, 3};
//+
Extrude {0, 0, -4} {
  Curve{2}; Curve{1}; Curve{4}; Curve{3}; 
}
//+
Curve Loop(5) = {2, 3, 4, -1};
//+
Surface(5) = {5};
//+
Curve Loop(7) = {12, 11, -9, 7};
//+
Surface(6) = {7};
//+
Physical Surface(1) = {5, 4, 3, 2, 1, 6};
//+
Surface Loop(1) = {5, 2, 1, 6, 4, 3};
//+
Volume(1) = {1};
//+
Physical Volume(2) = {1};

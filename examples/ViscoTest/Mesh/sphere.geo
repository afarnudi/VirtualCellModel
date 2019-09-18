// Gmsh project created on Fri Aug  9 21:01:05 2019
SetFactory("OpenCASCADE");
//+
Point(1) = {20, 0, 0, 1.0};
Point(2) = {-20, 0, 0, 1.0};
Point(3) = {0, 20, 0, 1.0};
Point(4) = {0, -20, 0, 1.0};
Point(5) = {0, 0, 20, 1.0};
Point(6) = {0, 0, -20, 1.0};
Point(7) = {0, 0, 0, 1.0};
//+
Circle(1) = {3, 7, 5};
//+
Circle(2) = {3, 7, 6};
//+
Circle(3) = {5, 7, 4};
//+
Circle(4) = {6, 7, 4};
//+
Circle(5) = {1, 7, 5};
//+
Circle(6) = {1, 7, 6};
//+
Circle(7) = {1, 7, 3};
//+
Circle(8) = {1, 7, 4};
//+
Circle(9) = {2, 7, 6};
//+
Circle(10) = {2, 7, 3};
//+
Circle(11) = {2, 7, 5};
//+
Circle(12) = {2, 7, 4};
//+
Curve Loop(1) = {1, -5, 7};
//+
Surface(1) = {1};
//+
Curve Loop(3) = {1, -11, 10};
//+
Surface(2) = {3};
//+
Curve Loop(5) = {3, -12, 11};
//+
Surface(3) = {5};
//+
Curve Loop(7) = {10, 2, -9};
//+
Surface(4) = {7};
//+
Curve Loop(9) = {7, 2, -6};
//+
Surface(5) = {9};
//+
Curve Loop(11) = {8, -4, -6};
//+
Surface(6) = {11};
//+
Curve Loop(13) = {3, -8, 5};
//+
Surface(7) = {13};
//+
Curve Loop(15) = {12, -4, -9};
//+
Surface(8) = {15};
//+
Physical Surface(1) = {1, 5, 4, 2, 3, 7, 6, 8};
//+
Surface Loop(1) = {1, 2, 3, 7, 6, 8, 4, 5};
//+
Volume(1) = {1};
//+
Physical Volume(2) = {1};

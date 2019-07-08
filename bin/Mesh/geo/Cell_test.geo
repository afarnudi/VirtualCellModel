r=30;

Point(1) = {0, 0, 0};
Point(2) = {r, 0, 0};
Point(3) = {-r, 0, 0};
Point(4) = {0, r, 0};
Point(5) = {0, -r, 0};
Point(6) = {0, 0, r};
Point(7) = {0, 0, -r};
//+
Circle(1) = {6, 1, 3};
//+
Circle(2) = {3, 1, 4};
//+
Circle(3) = {4, 1, 6};
//+
Circle(4) = {6, 1, 2};
//+
Circle(5) = {2, 1, 4};
//+
Circle(6) = {4, 1, 7};
//+
Circle(7) = {7, 1, 2};
//+
Circle(8) = {7, 1, 3};
//+
Circle(9) = {3, 1, 5};
//+
Circle(10) = {5, 1, 6};
//+
Circle(11) = {5, 1, 7};
//+
Circle(12) = {2, 1, 5};
//+
Curve Loop(1) = {3, 1, 2};
//+
Surface(1) = {1};
//+
Curve Loop(2) = {3, 4, 5};
//+
Surface(2) = {2};
//+
Curve Loop(3) = {5, 6, 7};
//+
Surface(3) = {3};
//+
Curve Loop(4) = {6, 8, 2};
//+
Surface(4) = {4};
//+
Curve Loop(5) = {8, 9, 11};
//+
Surface(5) = {5};
//+
Curve Loop(6) = {11, 7, 12};
//+
Surface(6) = {6};
//+
Curve Loop(7) = {12, 10, 4};
//+
Surface(7) = {7};
//+
Curve Loop(8) = {10, 1, 9};
//+
Surface(8) = {8};
//+
Surface Loop(1) = {6, 5, 4, 3, 2, 1, 8, 7};
//+
Volume(1) = {1};
//+
Physical Surface(1) = {6, 4, 3, 5, 8, 7, 1, 2};
//+
Physical Volume(2) = {1};

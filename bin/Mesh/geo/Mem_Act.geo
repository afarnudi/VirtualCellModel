s=4;
r=10;

Point(1) = {0, 0, 0, s};
Point(2) = {r, 0, 0, s};
Point(3) = {-r, 0, 0, s};
Point(4) = {0, r, 0, s};
Point(5) = {0, -r, 0, s};
Point(6) = {0, 0, r, s};
Point(7) = {0, 0, -r, s};

//+
Circle(1) = {6, 1, 3};
//+
Circle(2) = {3, 1, 7};
//+
Circle(3) = {7, 1, 2};
//+
Circle(4) = {2, 1, 6};
//+
Circle(5) = {6, 1, 4};
//+
Circle(6) = {4, 1, 7};
//+
Circle(7) = {7, 1, 5};
//+
Circle(8) = {5, 1, 6};
//+
Circle(9) = {3, 1, 4};
//+
Circle(10) = {4, 1, 2};
//+
Circle(11) = {2, 1, 5};
//+
Circle(12) = {5, 1, 3};
//+
Curve Loop(1) = {1, 9, -5};
//+
Surface(1) = {1};
//+
Curve Loop(2) = {2, -6, -9};
//+
Surface(2) = {2};
//+
Curve Loop(3) = {10, -3, -6};
//+
Surface(3) = {3};
//+
Curve Loop(4) = {10, 4, 5};
//+
Surface(4) = {4};
//+
Curve Loop(5) = {4, -8, -11};
//+
Surface(5) = {5};
//+
Curve Loop(6) = {11, -7, 3};
//+
Surface(6) = {6};
//+
Curve Loop(7) = {7, 12, 2};
//+
Surface(7) = {7};
//+
Curve Loop(8) = {1, -12, 8};
//+
Surface(8) = {8};
//+
Surface Loop(1) = {7, 6, 5, 4, 3, 2, 1, 8};
//+
Volume(1) = {1};
//+
Physical Surface(1) = {2, 3, 4, 5, 6, 7, 8, 1};
//+
Physical Volume(2) = {1};

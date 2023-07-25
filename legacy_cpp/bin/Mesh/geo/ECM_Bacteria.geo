s=3.0;
r=10;
l=30;

Point(1) = {0, 0, l/2, s};
Point(2) = {r, 0, l/2, s};
Point(3) = {-r, 0, l/2, s};
Point(4) = {0, r, l/2, s};
Point(5) = {0, -r, l/2, s};
Point(6) = {0, 0, l/2+r, s};

Point(7) = {0, 0, -l/2, s};
Point(8) = {r, 0, -l/2, s};
Point(9) = {-r, 0, -l/2, s};
Point(10) = {0, r, -l/2, s};
Point(11) = {0, -r, -l/2, s};
Point(12) = {0, 0, -l/2-r, s};//+
Circle(1) = {3, 1, 4};
//+
Circle(2) = {4, 1, 2};
//+
Circle(3) = {2, 1, 5};
//+
Circle(4) = {5, 1, 3};
//+
Circle(5) = {3, 1, 6};
//+
Circle(6) = {6, 1, 4};
//+
Circle(7) = {6, 1, 2};
//+
Circle(8) = {6, 1, 5};
//+
Circle(9) = {10, 7, 9};
//+
Circle(10) = {9, 7, 11};
//+
Circle(11) = {11, 7, 8};
//+
Circle(12) = {8, 7, 10};
//+
Circle(13) = {10, 7, 12};
//+
Circle(14) = {12, 7, 9};
//+
Circle(15) = {12, 7, 11};
//+
Circle(16) = {12, 7, 8};
//+
Line(17) = {10, 4};
//+
Line(18) = {3, 9};
//+
Line(19) = {11, 5};
//+
Line(20) = {2, 8};

//+
Curve Loop(1) = {1, -17, 9, -18};
//+
Surface(1) = {1};
//+
Curve Loop(2) = {17, 2, 20, 12};
//+
Surface(2) = {2};
//+
Curve Loop(3) = {20, -11, 19, -3};
//+
Surface(3) = {3};
//+
Curve Loop(4) = {18, 10, 19, 4};
//+
Surface(4) = {4};
//+
Curve Loop(5) = {10, -15, 14};
//+
Surface(5) = {5};
//+
Curve Loop(6) = {14, -9, 13};
//+
Surface(6) = {6};
//+
Curve Loop(7) = {13, 16, 12};
//+
Surface(7) = {7};
//+
Curve Loop(8) = {16, -11, -15};
//+
Surface(8) = {8};
//+
Curve Loop(9) = {6, -1, 5};
//+
Surface(9) = {9};
//+
Curve Loop(10) = {5, 8, 4};
//+
Surface(10) = {10};
//+
Curve Loop(11) = {8, -3, -7};
//+
Surface(11) = {11};
//+
Curve Loop(12) = {7, -2, -6};
//+
Surface(12) = {12};
//+
Physical Surface(1) = {3, 4, 1, 2, 8, 6, 5, 7, 10, 9, 12, 11};

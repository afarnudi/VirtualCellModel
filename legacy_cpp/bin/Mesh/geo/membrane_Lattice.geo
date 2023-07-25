s=3.0;
r=10;

Point(1) = {0, 0, 0, s};
Point(2) = {r, 0, 0, s};
Point(3) = {-r, 0, 0, s};
Point(4) = {0, r, 0, s};
Point(5) = {0, -r, 0, s};
Point(6) = {0, 0, r, s};
Point(7) = {0, 0, -r, s};

//+
Circle(1) = {3, 1, 4};
//+
Circle(2) = {6, 1, 3};
//+
Circle(3) = {6, 1, 4};
//+
Circle(4) = {6, 1, 2};
//+
Circle(5) = {2, 1, 4};
//+
Circle(6) = {2, 1, 7};
//+
Circle(7) = {7, 1, 4};
//+
Circle(8) = {7, 1, 3};
//+
Circle(9) = {7, 1, 5};
//+
Circle(10) = {3, 1, 5};
//+
Circle(11) = {5, 1, 6};
//+
Circle(12) = {5, 1, 2};
//+
Line Loop(13) = {1, -3, 2};
//+
Ruled Surface(14) = {13};
//+
Line Loop(15) = {6, 7, -5};
//+
Ruled Surface(16) = {15};
//+
Line Loop(17) = {5, -3, 4};
//+
Ruled Surface(18) = {17};
//+
Line Loop(19) = {7, -1, -8};
//+
Ruled Surface(20) = {19};
//+
Line Loop(21) = {11, 2, 10};
//+
Ruled Surface(22) = {21};
//+
Line Loop(23) = {10, -9, 8};
//+
Ruled Surface(24) = {23};
//+
Line Loop(25) = {4, -12, 11};
//+
Ruled Surface(26) = {25};
//+
Line Loop(27) = {6, 9, 12};
//+
Ruled Surface(28) = {27};
//+
Surface Loop(29) = {20, 16, 28, 24, 22, 26, 18, 14};
//+
Volume(30) = {29};
//+
Physical Surface(31) = {20, 16, 18, 14, 22, 24, 26, 28};
//+
Physical Volume(32) = {30};

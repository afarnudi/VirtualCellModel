s=2;
r=4;
l=10;

Point(1) = {l, 0, 0, s};
Point(2) = {l, r, r, s};
Point(3) = {l, r, -r, s};
Point(4) = {l, 2*r, 0, s};
Point(5) = {l, r, 0, s};
Point(6) = {l+r, r, 0, s};


Point(7) = {-l, 0, 0, s};
Point(8) = {-l, r, r, s};
Point(9) = {-l, r, -r, s};
Point(10) = {-l, 2*r, 0, s};
Point(11) = {-l, r, 0, s};
Point(12) = {-l-r, r, 0, s};
//+
Line(1) = {9, 3};
//+
Line(2) = {4, 10};
//+
Line(3) = {8, 2};
//+
Line(4) = {1, 7};
//+
Circle(5) = {2, 5, 4};
//+
Circle(6) = {4, 5, 3};
//+
Circle(7) = {3, 5, 1};
//+
Circle(8) = {1, 5, 2};
//+
Circle(9) = {2, 5, 6};
//+
Circle(10) = {6, 5, 1};
//+
Circle(11) = {6, 5, 3};
//+
Circle(12) = {6, 5, 4};
//+
Circle(13) = {9, 11, 10};
//+
Circle(14) = {10, 11, 8};
//+
Circle(15) = {8, 11, 7};
//+
Circle(16) = {7, 11, 9};
//+
Circle(17) = {9, 11, 12};
//+
Circle(18) = {12, 11, 7};
//+
Circle(19) = {12, 11, 8};
//+
Circle(20) = {12, 11, 10};
//+
Line Loop(21) = {2, 14, 3, 5};
//+
Ruled Surface(22) = {21};
//+
Line Loop(23) = {3, -8, 4, -15};
//+
Ruled Surface(24) = {23};
//+
Line Loop(25) = {16, 1, 7, 4};
//+
Ruled Surface(26) = {25};
//+
Line Loop(27) = {6, -1, 13, -2};
//+
Ruled Surface(28) = {27};
//+
Line Loop(29) = {14, -19, 20};
//+
Ruled Surface(30) = {29};
//+
Line Loop(31) = {15, -18, 19};
//+
Ruled Surface(32) = {31};
//+
Line Loop(33) = {18, 16, 17};
//+
Ruled Surface(34) = {33};
//+
Line Loop(35) = {20, -13, 17};
//+
Ruled Surface(36) = {35};
//+
Line Loop(37) = {6, -11, 12};
//+
Ruled Surface(38) = {37};
//+
Line Loop(39) = {12, -5, 9};
//+
Ruled Surface(40) = {39};
//+
Line Loop(41) = {9, 10, 8};
//+
Ruled Surface(42) = {41};
//+
Line Loop(43) = {10, -7, -11};
//+
Ruled Surface(44) = {43};
//+
Physical Surface(45) = {28, 22, 24, 26, 38, 40, 42, 44, 36, 34, 30, 32};

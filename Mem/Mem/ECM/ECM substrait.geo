s=4;
d=50;

Point(1) = {d, 0, d, s};
Point(2) = {d, 0, -d, s};
Point(3) = {-d, 0, d, s};
Point(4) = {-d, 0, -d, s};
//+
Line(1) = {4, 3};
//+
Line(2) = {1, 3};
//+
Line(3) = {1, 2};
//+
Line(4) = {2, 4};
//+
Line Loop(5) = {3, 4, 1, -2};
//+
Plane Surface(6) = {5};
//+
Physical Surface(7) = {6};

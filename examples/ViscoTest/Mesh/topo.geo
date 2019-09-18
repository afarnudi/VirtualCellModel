s=3.0;
r=10;
d=1;
inter=4;

Point(1) = {0, d, 2*r, s};
Point(2) = {0, d, -2*r, s};
Point(3) = {inter, d, 2*r, s};
Point(4) = {inter, d, -2*r, s};

Point(5) = {0, 0, 2*r, s};
Point(6) = {0, 0, -2*r, s};
Point(7) = {-inter, 0, 2*r, s};
Point(8) = {-inter, 0, -2*r, s};

Point(9) = {2*inter, d, 2*r, s};
Point(10) = {2*inter, d, -2*r, s};
Point(11) = {3*inter, d, 2*r, s};
Point(12) = {3*inter, d, -2*r, s};

Point(13) = {2*inter, 0, 2*r, s};
Point(14) = {2*inter, 0, -2*r, s};
Point(15) = {inter, 0, 2*r, s};
Point(16) = {inter, 0, -2*r, s};

Point(17) = {-2*inter, d, 2*r, s};
Point(18) = {-2*inter, d, -2*r, s};
Point(19) = {-inter, d, 2*r, s};
Point(20) = {-inter, d, -2*r, s};

Point(21) = {-2*inter, 0, 2*r, s};
Point(22) = {-2*inter, 0, -2*r, s};
Point(23) = {-3*inter, 0, 2*r, s};
Point(24) = {-3*inter, 0, -2*r, s};//+
//+
Line(1) = {23, 24};
//+
Line(2) = {21, 22};
//+
Line(3) = {17, 18};
//+
Line(4) = {19, 20};
//+
Line(5) = {7, 8};
//+
Line(6) = {5, 6};
//+
Line(7) = {1, 2};
//+
Line(8) = {3, 4};
//+
Line(9) = {15, 16};
//+
Line(10) = {13, 14};
//+
Line(11) = {9, 10};
//+
Line(12) = {11, 12};
//+
Line(13) = {23, 21};
//+
Line(14) = {21, 17};
//+
Line(15) = {17, 19};
//+
Line(16) = {19, 7};
//+
Line(17) = {7, 5};
//+
Line(18) = {5, 1};
//+
Line(19) = {1, 3};
//+
Line(20) = {3, 15};
//+
Line(21) = {15, 13};
//+
Line(22) = {13, 9};
//+
Line(23) = {9, 11};
//+
Line(24) = {24, 22};
//+
Line(25) = {22, 18};
//+
Line(26) = {18, 20};
//+
Line(27) = {20, 8};
//+
Line(28) = {8, 6};
//+
Line(29) = {6, 2};
//+
Line(30) = {2, 4};
//+
Line(31) = {4, 16};
//+
Line(32) = {16, 14};
//+
Line(33) = {14, 10};
//+
Line(34) = {10, 12};
//+
Curve Loop(1) = {1, 24, -2, -13};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {2, 25, -3, -14};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {15, 4, -26, -3};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {4, 27, -5, -16};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {5, 28, -6, -17};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {6, 29, -7, -18};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {7, 30, -8, -19};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {20, 9, -31, -8};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {9, 32, -10, -21};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {10, 33, -11, -22};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {11, 34, -12, -23};
//+
Plane Surface(11) = {11};
//+
Physical Surface(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

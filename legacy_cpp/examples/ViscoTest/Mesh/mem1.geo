// Gmsh project created on Sat Jul 27 12:06:13 2019
SetFactory("OpenCASCADE");
//+
a = 20;
b = 10;
Point(1) = {0, b, a, 1.0};
Point(2) = {0, -b, a, 1.0};
Point(3) = {b, 0, a, 1.0};
Point(4) = {-b, 0, a, 1.0};
Point(5) = {0, 0, a+b, 1.0};
Point(6) = {0, b, -a, 1.0};
Point(7) = {0, -b, -a, 1.0};
Point(8) = {b, 0, -a, 1.0};
Point(9) = {-b, 0, -a, 1.0};
Point(10) = {0, 0, -a-b, 1.0};
Point(11) = {0, 0, a, 1.0};
Point(12) = {0, 0, -a, 1.0};
//+
 

//+
Circle(1) = {1, 11, 4};
Circle(2) = {1, 11, 3};
Circle(3) = {2, 11, 4};
Circle(4) = {2, 11, 3};
Circle(5) = {1, 11, 5};
Circle(6) = {2, 11, 5};
Circle(7) = {3, 11, 5};
Circle(8) = {4, 11, 5};
Circle(9) = {6, 12, 8};
Circle(10) = {6, 12, 9};
Circle(11) = {7, 12, 8};
Circle(12) = {7, 12, 9};
Circle(13) = {6, 12, 10};
Circle(14) = {7, 12, 10};
Circle(15) = {8, 12, 10};
Circle(16) = {9, 12, 10};
//+
Line(17) = {1, 6};
//+
Line(18) = {3, 8};
//+
Line(19) = {4, 9};
//+
Line(20) = {2, 7};
//+
Curve Loop(1) = {2, 7, -5};
//+
Surface(1) = {1};
//+
Curve Loop(3) = {1, 8, -5};
//+
Surface(2) = {3};
//+
Curve Loop(5) = {6, -8, -3};
//+
Surface(3) = {5};
//+
Curve Loop(7) = {7, -6, 4};
//+
Surface(4) = {7};
//+
Curve Loop(9) = {2, 18, -9, -17};
//+
Surface(5) = {9};
//+
Curve Loop(11) = {4, 18, -11, -20};
//+
Surface(6) = {11};
//+
Curve Loop(13) = {20, 12, -19, -3};
//+
Surface(7) = {13};
//+
Curve Loop(15) = {19, -10, -17, 1};
//+
Surface(8) = {15};
//+
Curve Loop(17) = {16, -13, 10};
//+
Surface(9) = {17};
//+
Curve Loop(19) = {16, -14, 12};
//+
Surface(10) = {19};
//+
Curve Loop(21) = {14, -15, -11};
//+
Surface(11) = {21};
//+
Curve Loop(23) = {13, -15, -9};
//+
Surface(12) = {23};

//+
Physical Surface(1) = {1, 2, 3, 4, 6, 5, 8, 7, 10, 11, 12, 9};
//+
Surface Loop(1) = {7, 6, 4, 1, 2, 3, 8, 9, 10, 11, 12, 5};
//+
Volume(1) = {1};
//+
Physical Volume(2) = {1};

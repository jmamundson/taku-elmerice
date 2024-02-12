L=10.0;
reslt=L/10.0;
//+
Point(1) = {-L/2, 0, 0, reslt};
//+
Point(2) = {L/2, -0, 0, reslt};
//+
Line(1) = {1, 2};
//+
Extrude {0, 1, 0} {
  Curve{1}; Layers{10}; Recombine;
}
//+
Physical Surface(6) = {5};
//+
Physical Curve(7) = {3};
//+
Physical Curve(8) = {4};
//+
Physical Curve(9) = {1};
//+
Physical Curve(10) = {2};

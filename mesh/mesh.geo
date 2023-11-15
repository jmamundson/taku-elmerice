DefineConstant[ frontres = {100, Name "frontres"}]; // mesh size near glacier terminus
DefineConstant[ backres = { 500, Name "backres"}]; // mesh size near ice divide
Point(1) = {0, 0, 0, backres}; // coordinate of bedrock at ice divide
Point(2) = {40000, 0, 0, frontres}; // coordinate of bedrock at glacier terminus
Line(1) = {1, 2}; // create line from points 1 and 2
Extrude {0, 1000, 0} {
  Line{1}; Layers{10}; Recombine;
}
Physical Surface(6) = {5};
Physical Line(7) = {1};
Physical Line(8) = {2};
Physical Line(9) = {3};
Physical Line(10) = {4};

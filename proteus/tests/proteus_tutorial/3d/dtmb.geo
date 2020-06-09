//+
SetFactory("OpenCASCADE");
Geometry.Tolerance = 1e-2;
Geometry.OCCSewFaces = 1;
Geometry.OCCFixSmallEdges = 1;
Geometry.OCCFixSmallFaces = 1;
Geometry.OCCMakeSolids = 1;
//+
Merge "dtmb.igs";
//+
bv = newv;
//+
Box(bv) = {-2, -4, -1, 20, 8, 2};
//+
tv = newv;
//+
BooleanDifference(tv) = {Volume{bv}; Delete; }{ Volume{1}; Delete; };
//+
Physical Surface(7) = {14, 11, 9, 7, 13, 10, 8, 12, 15, 18, 17, 19, 16};
//+
Physical Volume(1) = {tv};
//+
Physical Point(1) = {1, 3, 7, 6};
//+
Physical Point(2) = {2, 4, 8, 5};
//+
Physical Curve(2) = {4, 10, 12, 5};
//+
Physical Curve(1) = {2, 8, 9, 7};
//+
Physical Surface(1) = {3};
//+
Physical Surface(2) = {5};
//+
Physical Surface(3) = {1};
//+
Physical Surface(4) = {6};
//+
Physical Surface(5) = {2};
//+
Physical Surface(6) = {4};
//+
Physical Curve(3) = {1, 3};
//+
Physical Curve(4) = {6, 11};

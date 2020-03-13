
// Points
Point(1) = {0,0,0};
Point(2) = {1,0,0};
Point(3) = {1,1,0};
Point(4) = {0,1,0};
Point(5) = {0,0,1};
Point(6) = {1,0,1};
Point(7) = {1,1,1};
Point(8) = {0,1,1};
Point(9) = {0.4,0.4,0.41};
Point(10) = {0.4,0.6,0.41};
Point(11) = {0.6,0.6,0.41};
Point(12) = {0.6,0.4,0.41};
Point(13) = {0.4,0.4,0.61};
Point(14) = {0.4,0.6,0.61};
Point(15) = {0.6,0.6,0.61};
Point(16) = {0.6,0.4,0.61};

// Lines

// Surfaces
Line(1) = {4,1};
Line(2) = {1,2};
Line(3) = {2,3};
Line(4) = {3,4};
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Line(6) = {8,5};
Line(7) = {5,6};
Line(8) = {6,7};
Line(9) = {7,8};
Line Loop(2) = {6, 7, 8, 9};
Plane Surface(2) = {2};
Line(11) = {5,1};
Line(12) = {1,2};
Line(13) = {2,6};
Line(14) = {6,5};
Line Loop(3) = {11, 12, 13, 14};
Plane Surface(3) = {3};
Line(16) = {6,2};
Line(17) = {2,3};
Line(18) = {3,7};
Line(19) = {7,6};
Line Loop(4) = {16, 17, 18, 19};
Plane Surface(4) = {4};
Line(21) = {8,4};
Line(22) = {4,3};
Line(23) = {3,7};
Line(24) = {7,8};
Line Loop(5) = {21, 22, 23, 24};
Plane Surface(5) = {5};
Line(26) = {5,1};
Line(27) = {1,4};
Line(28) = {4,8};
Line(29) = {8,5};
Line Loop(6) = {26, 27, 28, 29};
Plane Surface(6) = {6};
Line(31) = {12,9};
Line(32) = {9,10};
Line(33) = {10,11};
Line(34) = {11,12};
Line Loop(7) = {31, 32, 33, 34};
Plane Surface(7) = {7};
Line(36) = {14,10};
Line(37) = {10,11};
Line(38) = {11,15};
Line(39) = {15,14};
Line Loop(8) = {36, 37, 38, 39};
Plane Surface(8) = {8};
Line(41) = {15,11};
Line(42) = {11,12};
Line(43) = {12,16};
Line(44) = {16,15};
Line Loop(9) = {41, 42, 43, 44};
Plane Surface(9) = {9};
Line(46) = {16,12};
Line(47) = {12,9};
Line(48) = {9,13};
Line(49) = {13,16};
Line Loop(10) = {46, 47, 48, 49};
Plane Surface(10) = {10};
Line(51) = {13,9};
Line(52) = {9,10};
Line(53) = {10,14};
Line(54) = {14,13};
Line Loop(11) = {51, 52, 53, 54};
Plane Surface(11) = {11};
Line(56) = {16,13};
Line(57) = {13,14};
Line(58) = {14,15};
Line(59) = {15,16};
Line Loop(12) = {56, 57, 58, 59};
Plane Surface(12) = {12};

// Volumes
Surface Loop(13) = {1, 2, 3, 4, 5, 6};
Surface Loop(14) = {7, 8, 9, 10, 11, 12};
Volume(1) = {13, 14};

// Physical Groups
Physical Point(1) = {1, 2, 3, 4};
Physical Point(6) = {5, 6, 7, 8};
Physical Point(9) = {9, 10, 11, 12};
Physical Point(14) = {13, 14, 15, 16};
Physical Surface(1) = {1};
Physical Surface(6) = {2};
Physical Surface(5) = {3};
Physical Surface(4) = {4};
Physical Surface(3) = {5};
Physical Surface(2) = {6};
Physical Surface(9) = {7};
Physical Surface(10) = {8};
Physical Surface(11) = {9};
Physical Surface(12) = {10};
Physical Surface(13) = {11};
Physical Surface(14) = {12};
Physical Volume(1) = {1};

// Other Options
Coherence;

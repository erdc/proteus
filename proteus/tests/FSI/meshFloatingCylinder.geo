
// Points
Point(1) = {0,0,0};
Point(2) = {1,0,0};
Point(3) = {1,1,0};
Point(4) = {0,1,0};
Point(5) = {0.4,0.41,0};
Point(6) = {0.6,0.41,0};
Point(7) = {0.6,0.61,0};
Point(8) = {0.4,0.61,0};

// Lines
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,5};

// Surfaces
Line Loop(1) = {4, 1, 2, 3};
Line Loop(2) = {8, 5, 6, 7};
Plane Surface(1) = {1, 2};

// Volumes

// Physical Groups
Physical Point(1) = {1, 2};
Physical Point(3) = {3, 4};
Physical Point(6) = {5, 6};
Physical Point(8) = {7, 8};
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
Physical Line(6) = {5};
Physical Line(7) = {6};
Physical Line(8) = {7};
Physical Line(9) = {8};
Physical Surface(1) = {1};

// Other Options
Coherence;

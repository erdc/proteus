/* geometric parameters */
L = 1;
ow = 0.33;
oh = 0.2;
och = 0.25;
/* meshing parameter */
size = 0.1;
/* corner points of the tank */
Point(0) = {0,0,0,size};
Point(1) = {L,0,0,size};
Point(2) = {L,L,0,size};
Point(3) = {0,L,0,size};
Point(4) = {0,0,L,size};
Point(5) = {L,0,L,size};
Point(6) = {L,L,L,size};
Point(7) = {0,L,L,size};
/* edges of the tank */
Line( 1) = {0,1};
Line( 2) = {1,2};
Line( 3) = {2,3};
Line( 4) = {3,0};
Line( 5) = {4,5};
Line( 6) = {5,6};
Line( 7) = {6,7};
Line( 8) = {7,4};
Line( 9) = {0,4};
Line(10) = {1,5};
Line(11) = {2,6};
Line(12) = {3,7};
/* faces of the tank */
/* the four "side" faces */
Line Loop(1) = {  1, 10,- 5,- 9};
Plane Surface(1) = {1};
Line Loop(2) = {  2, 11,- 6,-10};
Plane Surface(2) = {2};
Line Loop(3) = {  3, 12,- 7,-11};
Plane Surface(3) = {3};
Line Loop(4) = {  4,  9,- 8,-12};
Plane Surface(4) = {4};
/* the two "cap" faces */
Line Loop(5) = {- 4,- 3,- 2,- 1};
Plane Surface(5) = {5};
Line Loop(6) = {  5,  6,  7,  8};
Plane Surface(6) = {6};
/* the shell for the tank's outer boundary */
Surface Loop(1) = { 1, 2, 3, 4, 5, 6};
/* corner points of the object */
Point( 8) = {L/2-ow/2,L/2-ow/2,och-oh/2,size};
Point( 9) = {L/2+ow/2,L/2-ow/2,och-oh/2,size};
Point(10) = {L/2+ow/2,L/2+ow/2,och-oh/2,size};
Point(11) = {L/2-ow/2,L/2+ow/2,och-oh/2,size};
Point(12) = {L/2-ow/2,L/2-ow/2,och+oh/2,size};
Point(13) = {L/2+ow/2,L/2-ow/2,och+oh/2,size};
Point(14) = {L/2+ow/2,L/2+ow/2,och+oh/2,size};
Point(15) = {L/2-ow/2,L/2+ow/2,och+oh/2,size};
/* edges of the object */
Line(12 +  1) = {8 + 0,8 + 1};
Line(12 +  2) = {8 + 1,8 + 2};
Line(12 +  3) = {8 + 2,8 + 3};
Line(12 +  4) = {8 + 3,8 + 0};
Line(12 +  5) = {8 + 4,8 + 5};
Line(12 +  6) = {8 + 5,8 + 6};
Line(12 +  7) = {8 + 6,8 + 7};
Line(12 +  8) = {8 + 7,8 + 4};
Line(12 +  9) = {8 + 0,8 + 4};
Line(12 + 10) = {8 + 1,8 + 5};
Line(12 + 11) = {8 + 2,8 + 6};
Line(12 + 12) = {8 + 3,8 + 7};
/* faces of the object */
/* the four "side" faces */
Line Loop(6 + 1) = {
 (12 +  1),
 (12 + 10),
-(12 +  5),
-(12 +  9)};
Plane Surface(6 + 1) = {6 + 1};
Line Loop(6 + 2) = {
 (12 +  2),
 (12 + 11),
-(12 +  6),
-(12 + 10)};
Plane Surface(6 + 2) = {6 + 2};
Line Loop(6 + 3) = {
 (12 +  3),
 (12 + 12),
-(12 +  7),
-(12 + 11)};
Plane Surface(6 + 3) = {6 + 3};
Line Loop(6 + 4) = {
 (12 +  4),
 (12 +  9),
-(12 +  8),
-(12 + 12)};
Plane Surface(6 + 4) = {6 + 4};
/* the two "cap" faces */
Line Loop(6 + 5) = {
-(12 +  4),
-(12 +  3),
-(12 +  2),
-(12 +  1)};
Plane Surface(6 + 5) = {6 + 5};
Line Loop(6 + 6) = {
 (12 +  5),
 (12 +  6),
 (12 +  7),
 (12 +  8)};
Plane Surface(6 + 6) = {6 + 6};
/* the shell for the objects's outer boundary */
Surface Loop(2) = {6+1,6+2,6+3,6+4,6+5,6+6};
/* the actual tank minus the object in it */
Volume(1) = {  1, -2};

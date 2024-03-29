Point(1) = { 8.0, 0.0, 8.0, 1.0 };
Point(2) = { 8.0, 8.0, 8.0, 1.0 };
Point(3) = { 0.0, 8.0, 8.0, 1.0 };
Point(4) = { 0.0, 0.0, 8.0, 1.0 };
Point(5) = { 8.0, 8.0, 0.0, 1.0 };
Point(6) = { 8.0, 0.0, 0.0, 1.0 };
Point(7) = { 0.0, 0.0, 0.0, 1.0 };
Point(8) = { 0.0, 8.0, 0.0, 1.0 };
Line(1) = { 1, 2 };
Line(2) = { 2, 3 };
Line(3) = { 3, 4 };
Line(4) = { 1, 4 };
Line(5) = { 5, 6 };
Line(6) = { 6, 7 };
Line(7) = { 7, 8 };
Line(8) = { 5, 8 };
Line(9) = { 2, 5 };
Line(10) = { 1, 6 };
Line(11) = { 4, 7 };
Line(12) = { 3, 8 };
Line Loop(1) = { 1, 2, 3, -4 };
Line Loop(2) = { 5, 6, 7, -8 };
Line Loop(3) = { -5, -9, -1, 10 };
Line Loop(4) = { -7, -11, -3, 12 };
Line Loop(5) = { 8, -12, -2, 9 };
Line Loop(6) = { -6, -10, 4, 11 };
Plane Surface(1) = { 1 };
Plane Surface(2) = { 2 };
Plane Surface(3) = { 3 };
Plane Surface(4) = { 4 };
Plane Surface(5) = { 5 };
Plane Surface(6) = { 6 };
Surface Loop(1) = { 1, 2, 3, 4, 5, 6 };
Physical Surface("Ceiling$0.16,0.16,0.16,0.16,0.16") = { 1 };
Physical Surface("Floor$0.16,0.16,0.16,0.16,0.16") = { 2 };
Physical Surface("WallRight$0.16,0.16,0.16,0.16,0.16") = { 3 };
Physical Surface("WallLeft$0.16,0.16,0.16,0.16,0.16") = { 4 };
Physical Surface("WallBack$0.16,0.16,0.16,0.16,0.16") = { 5 };
Physical Surface("Wallfront$0.16,0.16,0.16,0.16,0.16") = { 6 };
Volume( 1 ) = { 1 };
Physical Volume("region_1") = { 1 };
Physical Line ("default") = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
Mesh.Algorithm = 6;
Mesh.Algorithm3D = 1; // Delaunay3D, works for boundary layer insertion.
Mesh.Optimize = 1; // Gmsh smoother, works with boundary layers (netgen version does not).
Mesh.CharacteristicLengthFromPoints = 1;
// Recombine Surface "*";


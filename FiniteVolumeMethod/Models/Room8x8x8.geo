Point(1) = { 0.0, 0.0, 0.0, 1.0 };
Point(2) = { 8.0, 0.0, 0.0, 1.0 };
Point(3) = { 8.0, 0.0, 8.0, 1.0 };
Point(4) = { 0.0, 0.0, 8.0, 1.0 };
Point(5) = { 8.0, 8.0, 0.0, 1.0 };
Point(6) = { 8.0, 8.0, 8.0, 1.0 };
Point(7) = { 0.0, 8.0, 8.0, 1.0 };
Point(8) = { 0.0, 8.0, 0.0, 1.0 };
Line(1) = { 1, 2 };
Line(2) = { 2, 3 };
Line(3) = { 3, 4 };
Line(4) = { 1, 4 };
Line(5) = { 2, 5 };
Line(6) = { 5, 6 };
Line(7) = { 3, 6 };
Line(8) = { 4, 7 };
Line(9) = { 6, 7 };
Line(10) = { 1, 8 };
Line(11) = { 7, 8 };
Line(12) = { 5, 8 };
Line Loop(1) = { 1, 2, 3, -4 };
Line Loop(2) = { 5, 6, -7, -2 };
Line Loop(3) = { -8, -3, 7, 9 };
Line Loop(4) = { -10, 4, 8, 11 };
Line Loop(5) = { 10, -12, -5, -1 };
Line Loop(6) = { 12, -11, -9, -6 };
Plane Surface(1) = { 1 };
Plane Surface(2) = { 2 };
Plane Surface(3) = { 3 };
Plane Surface(4) = { 4 };
Plane Surface(5) = { 5 };
Plane Surface(6) = { 6 };
Surface Loop(1) = { 1, 2, 3, 4, 5, 6 };
Physical Surface("Sides") = { 1, 2, 4, 6 };
Physical Surface("TB") = { 3, 5 };
Volume( 1 ) = { 1 };
Physical Volume("Room8x8x8") = { 1 };
Physical Line ("default") = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
Mesh.Algorithm = 6;
Mesh.Algorithm3D = 1; // Delaunay3D, works for boundary layer insertion.
Mesh.Optimize = 1; // Gmsh smoother, works with boundary layers (netgen version does not).
Mesh.CharacteristicLengthFromPoints = 1;
// Recombine Surface "*";
Mesh.RemeshAlgorithm = 1; // automatic

//+
SetFactory("Built-in");

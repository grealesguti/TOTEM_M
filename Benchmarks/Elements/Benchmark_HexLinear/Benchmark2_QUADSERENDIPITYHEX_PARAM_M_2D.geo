// Gmsh project created on Mon Nov 06 12:14:39 2023
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, 0.0014, 0.001524, 0};
//+
Physical Surface("Surface", 5) = {1};
//+
Physical Curve("Ymin", 6) = {1};
//+
Physical Curve("Ymax", 7) = {3};
//+
Physical Curve("Xmin", 8) = {4};
//+
Physical Curve("Xmax", 9) = {2};
//+
Transfinite Surface {1};

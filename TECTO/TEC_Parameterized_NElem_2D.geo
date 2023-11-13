// Gmsh project created on Mon Nov 06 18:40:42 2023
Copper_Height=0.1;
Ceramic_Height=0.2;
Semiconductor_Height=1.2;
Semiconductor_Width=1;
Copper_bridge=0.2;
electrode_width=0.1;

element_copper=5;
element_ceramic=5;
element_semi_h=40;
element_semi_w=40;
element_bridge=5;


SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, Semiconductor_Width, Semiconductor_Height, 0};
//+
Extrude {0, Copper_Height, 0} {
  Curve{3}; 
}
//+
Extrude {Copper_bridge, 0, 0} {
  Curve{5}; 
}
//+
Extrude {Semiconductor_Width, 0, 0} {
  Curve{10}; 
}
//+
Extrude {0, -Semiconductor_Height, 0} {
  Curve{11}; 
}
//+
Extrude {0, -Copper_Height, 0} {
  Curve{16}; Curve{1}; 
}
//+
Extrude {electrode_width, 0, 0} {
  Curve{18}; 
}
//+
Extrude {-electrode_width, 0, 0} {
  Curve{20}; 
}
//+
Extrude {0, -Ceramic_Height, 0} {
  Curve{22}; Curve{27}; Curve{19}; Curve{24}; 
}
//+
Extrude {0, Ceramic_Height, 0} {
  Curve{7}; Curve{9}; Curve{12}; 
}
//+
Transfinite Curve {2, 4, 14, 14, 15, 15} = element_semi_h+1 Using Progression 1;
//+
Transfinite Curve {20, 28, 21, 17, 18, 25, 13, 10, 5, 6, 6} = element_copper+1 Using Progression 1;
//+
Transfinite Curve {32, 29, 30, 34, 34, 35, 35, 37, 37, 44, 44, 42, 42, 39, 40} = element_ceramic+1 Using Progression 1;
//+
Transfinite Curve {41, 7, 3, 1, 22, 31, 36, 19, 19, 16, 11, 12, 12, 45, 45} = element_semi_w+1 Using Progression 1;
//+
Transfinite Curve {43, 9, 9, 8, 8} = element_bridge+1 Using Progression 1;
//+
Transfinite Surface {9};
//+
Transfinite Surface {11};
//+
Transfinite Surface {7};
//+
Transfinite Surface {10};
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {14};
//+
Transfinite Surface {15};
//+
Transfinite Surface {3};
//+
Transfinite Surface {4};
//+
Transfinite Surface {16};
//+
Transfinite Surface {5};
//+
Transfinite Surface {6};
//+
Transfinite Surface {12};
//+
Transfinite Surface {13};
//+
Transfinite Surface {8};
//+
Physical Surface("Copper", 46) = {9, 7, 6, 8, 4, 3, 2};
//+
Physical Surface("SemiconductorL", 47) = {1};
//+
Physical Surface("SemiconductorR", 48) = {5};
//+
Physical Surface("Ceramic", 49) = {14, 15, 16, 12, 13, 10, 11};
//+
Physical Surface("Surface", 50) = {11, 10, 7, 9, 1, 2, 14, 15, 3, 4, 16, 5, 6, 12, 13, 8};
//+
Physical Curve("Ceramic_Top", 51) = {41, 43, 45};
//+
Physical Curve("Ceramic_Bottom", 52) = {31, 33, 36, 38};
//+
Physical Curve("Electrode_MinX", 53) = {28};
//+
Physical Curve("Electrode_MaxX", 54) = {25};
//+
Physical Curve("Contact_Ceramic", 55) = {27, 22, 19, 24, 12, 7, 9};
//+
Physical Surface("TO", 56) = {5, 1};

//+
Physical Curve("Lateral_Force_L", 57) = {40};

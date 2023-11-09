// Gmsh project created on Mon Nov 06 23:03:43 2023
Copper_Height=0.1;
Solder_Height=0.05;
Semiconductor_Height=1.2;
Semiconductor_separation=0.1;
Semiconductor_Block=Semiconductor_Height-2*Semiconductor_separation;
Semiconductor_width=1;
Ceramic_Height=0.2;
Copper_bridge=0.2;

SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, Semiconductor_width, Semiconductor_Block, 0};
//+
Extrude {0, Semiconductor_separation, 0} {
  Curve{3}; 
}
//+
Extrude {0, Solder_Height, 0} {
  Curve{7}; 
}
//+
Extrude {0, Copper_Height, 0} {
  Curve{10}; 
}
//+
Extrude {Copper_bridge, 0, 0} {
  Curve{11}; 
}
//+
Extrude {Semiconductor_width, 0, 0} {
  Curve{16}; 
}
//+
Extrude {0, -Solder_Height, 0} {
  Curve{17}; 
}
//+
Extrude {0, -Semiconductor_separation, 0} {
  Curve{22}; 
}
//+
Extrude {0, -Semiconductor_Block, 0} {
  Curve{25}; 
}
//+
Extrude {0, -Semiconductor_separation, 0} {
  Curve{28}; Curve{1}; 
}
//+
Extrude {0, -Solder_Height, 0} {
  Curve{34}; Curve{31}; 
}
//+
Extrude {0, -Copper_Height, 0} {
  Curve{37}; Curve{40}; 
}
//+
Extrude {-0.2, 0, 0} {
  Curve{41}; 
}
//+
Extrude {0.2, 0, 0} {
  Curve{45}; 
}
//+
Extrude {0, -Ceramic_Height, 0} {
  Curve{48}; Curve{43}; Curve{46}; Curve{51}; 
}
//+
Extrude {0, Ceramic_Height, 0} {
  Curve{13}; Curve{15}; Curve{18}; 
}//+
Extrude {-0.2, 0, 0} {
  Curve{4}; 
}
//+
Extrude {+0.2, 0, 0} {
  Curve{27}; 
}

Copper_elements=3;
Solder_elements=2;
Copper_bridge_elements=3;
Semiconductor_block_elements_height=20;
Semiconductor_block_elements_width=20;

Ceramic_elements=5;
//+
Transfinite Curve {72, 4, 2, 26, 27, 75} = Semiconductor_block_elements_height+1 Using Progression 1;
//+
Transfinite Curve {65, 13, 10, 7, 7, 3, 1, 34, 37, 37, 43, 43, 57, 57, 60, 60, 46, 40, 40, 31, 31, 28, 25, 25, 22, 22, 17, 18, 18, 69} = Semiconductor_block_elements_width+1 Using Progression 1;
//+
Transfinite Curve {49, 41, 42, 44, 45, 52, 19, 16, 11, 12} = Copper_elements+1 Using Progression 1;
//+
Transfinite Curve {64, 63, 66, 68, 68, 61, 59, 58, 56, 53, 54} = Ceramic_elements+1 Using Progression 1;
//+
Transfinite Curve {35, 36, 38, 39, 39, 21, 20, 8, 9} = Solder_elements+1 Using Progression 1;
//+
//+
Transfinite Surface {18};
//+
Transfinite Surface {16};
//+
Transfinite Surface {14};
//+
Transfinite Surface {19};
//+
Transfinite Surface {12};
//+
Transfinite Surface {11};
//+
Transfinite Surface {1};
//+
Transfinite Surface {25};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3};
//+
Transfinite Surface {4};
//+
Transfinite Surface {22};
//+
Transfinite Surface {23};
//+
Transfinite Surface {5};
//+
Transfinite Surface {24};
//+
Transfinite Surface {6};
//+
Transfinite Surface {7};
//+
Transfinite Surface {9};
//+
Transfinite Surface {8};
//+
Transfinite Surface {26};
//+
Transfinite Surface {10};
//+
Transfinite Surface {13};
//+
Transfinite Surface {15};
//+
Transfinite Surface {20};
//+
Transfinite Surface {21};
//+
Transfinite Surface {17};
//+
Physical Surface("SemiconductorL", 76) = {1, 11, 25, 2};
//+
Physical Surface("SemiconductorR", 77) = {9, 26, 8, 10};
//+
Physical Surface("Solder", 78) = {12, 13, 7, 3};
//+
Physical Surface("Copper", 79) = {4, 5, 6, 15, 17, 14, 16};
//+
Physical Surface("Ceramic", 80) = {18, 19, 20, 21, 24, 23, 22};
//+
Physical Surface("Surface", 81) = {18, 19, 14, 16, 12, 11, 1, 25, 2, 3, 4, 22, 23, 5, 24, 6, 7, 8, 9, 26, 10, 13, 15, 20, 21, 17};
//+
Physical Surface("TO", 82) = {2, 1, 25, 11, 10, 9, 26, 8};
//+
Physical Curve("Ceramic_Bottom", 83) = {55, 57, 60, 62};
//+
Physical Curve("Ceramic_Top", 84) = {65, 67, 69};
//+
Physical Curve("Electrode_MaxX", 85) = {52};
//+
Physical Curve("Electrode_MinX", 86) = {49};
//+
Physical Curve("Contact_Ceramic", 87) = {13, 15, 18, 46, 51, 43, 48};
//+
Transfinite Curve {70, 71, 74, 73} = 5 Using Progression 1;
//+
Transfinite Curve {67, 15, 15, 14, 14} = 5 Using Progression 1;
//+
Transfinite Curve {47, 48, 55, 50, 51, 62} = 5 Using Progression 1;

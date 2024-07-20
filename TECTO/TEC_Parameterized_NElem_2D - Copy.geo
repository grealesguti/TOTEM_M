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
  Curve{3}; Recombine;
}
//+
Extrude {Copper_bridge, 0, 0} {
  Curve{5}; Recombine;
}
//+
Extrude {Semiconductor_Width, 0, 0} {
  Curve{10}; Recombine;
}
//+
Extrude {0, -Semiconductor_Height, 0} {
  Curve{11}; Recombine;
}
//+
Extrude {0, -Copper_Height, 0} {
  Curve{16}; Curve{1}; 
}
//+
Extrude {electrode_width, 0, 0} {
  Curve{18}; Recombine;
}
//+
Extrude {-electrode_width, 0, 0} {
  Curve{20}; Recombine;
}
//+
Extrude {0, -Ceramic_Height, 0} {
  Curve{22}; Curve{27}; Curve{19}; Curve{24}; Recombine;
}
//+
Extrude {0, Ceramic_Height, 0} {
  Curve{7}; Curve{9}; Curve{12}; Recombine;
}
Extrude {+Ceramic_Height, 0, 0} {
  Curve{39};  Recombine;
}
//+
Extrude {0, -Copper_Height, 0} {
  Curve{23};  Recombine;
}
//+
Extrude {0, Ceramic_Height, 0} {
  Curve{33};  Recombine;
}
Extrude {0, -0.2, 0} {
  Curve{45}; Recombine;
}
Recursive Delete {
  Surface{15};Surface{16}; Surface{8}; Surface{11}; 
}//+
Coherence;
//+
Transfinite Curve {6, 5, 44, 46, 67, 48, 47, 51, 55, 50} = element_copper+1 Using Progression 1;
//+
Transfinite Curve {71, 66, 39, 40, 69, 56, 57, 59, 59, 60, 62} = element_ceramic+1 Using Progression 1;
//+
Transfinite Curve {41, 7, 3, 3, 1, 52, 58, 61, 49, 16, 11, 11, 45, 72} = element_semi_w+1 Using Progression 1;
//+
Transfinite Curve {4, 2, 14, 15, 15} = element_semi_h+1 Using Progression 1;
//+
Transfinite Curve {65, 43, 43, 42, 42} = element_bridge+1 Using Progression 1;
//+
Transfinite Curve {68, 63, 64, 53, 54, 70, 70} = 4 Using Progression 1;
//+
Transfinite Surface {1};
//+
Transfinite Surface {5};
//+
Transfinite Surface {4};
//+
Transfinite Surface {20};
//+
Transfinite Surface {17};
//+
Transfinite Surface {3};
//+
Transfinite Surface {2};
//+
Transfinite Surface {14};
//+
Transfinite Surface {7};
//+
Transfinite Surface {10};
//+
Transfinite Surface {6};
//+
Transfinite Surface {12};
//+
Transfinite Surface {13};
//+
Transfinite Surface {18};
//+
Transfinite Surface {9};
//+
Transfinite Surface {19};


//+
Physical Surface("Copper", 73) = {9, 7, 6, 18, 4, 3, 2};
//+
Physical Surface("SemiconductorP+", 74) = {1};
//+
Physical Surface("SemiconductorN-", 75) = {5};
//+
Physical Surface("Ceramic", 76) = {14, 17, 20, 12, 13, 10, 19};

//+
Physical Curve("Ceramic_Top", 77) = {41, 65, 72};
//+
Physical Curve("Ceramic_Bottom", 78) = {58, 70, 61, 64};
//+
Physical Curve("ElectrodeMinX", 79) = {55};
//+
Physical Curve("ElectrodeMaxX", 80) = {67};


//+
Physical Surface("TO", 82) = {1, 5};
//+
Physical Surface("Volume", 83) = {1, 2, 14, 17, 20, 4, 3, 5, 6, 18, 13, 12, 10, 7, 9, 19};
//+
Physical Curve("Contact_CeramicCopper", 84) = {7, 43, 45, 49, 63, 52, 54};
//+

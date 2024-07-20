// Gmsh project created on Mon Nov 06 18:40:42 2023
Copper_Height=0.1;
Ceramic_Height=0.2;
Semiconductor_Height=1.2;
Semiconductor_Width=1;
Copper_bridge=0.2;
electrode_width=0.1;
solder_height=0.05;


element_copper=5;
element_ceramic=5;
element_semi_h=40;
element_semi_w=40;
element_bridge=5;
element_solder=2;



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
}//+//+
Translate {0, -solder_height, 0} {
  Surface{10}; Surface{7}; Surface{9}; Surface{19}; Surface{6}; Surface{12}; Surface{13}; Surface{18}; 
}
//+
Translate {0, solder_height, 0} {
  Surface{2}; Surface{14}; Surface{17}; Surface{3}; Surface{20}; Surface{4}; 
}
//+
Extrude {0, -solder_height, 0} {
  Curve{88};  Recombine;
}
//+
Extrude {0, -solder_height, 0} {
  Curve{103};  Recombine;
}
//+//+
Extrude {0, +solder_height, 0} {
  Curve{74};  Recombine;
}
//+
Extrude {0, +solder_height, 0} {
  Curve{64};  Recombine;
}
Coherence;
Extrude {0, -solder_height, 0} {
  Curve{104};  Recombine;
}
//+
Extrude {0, -solder_height, 0} {
  Curve{116};  Recombine;
}
//+
Extrude {0, +solder_height, 0} {
  Curve{106};  Recombine;
}
Recursive Delete {
  Surface{24};Surface{23}; Surface{21}; 
}//+

Extrude {-0.2, 0, 0} {
  Curve{129}; Recombine;
}
Coherence;
//+//+
Transfinite Curve {107, 105, 113, 114} = element_semi_h+1 Using Progression 1;
//+
Transfinite Curve {169, 168, 145, 146, 167, 166, 164, 165} = element_solder+1 Using Progression 1;
//+
Transfinite Curve {159, 154, 155, 151, 152, 162, 112, 110, 147, 148} = element_copper+1 Using Progression 1;
//+
Transfinite Curve {139, 126, 160, 129, 130, 132, 141, 136, 90, 91} = element_ceramic+1 Using Progression 1;
//+
Transfinite Curve {92, 89, 149, 106, 104, 156, 156, 157, 157, 161, 161, 131, 131, 119, 119, 153, 153, 116, 116, 115, 115, 103, 103, 111, 142, 142} = element_semi_w+1 Using Progression 1;
//+
Transfinite Curve {170, 171, 150, 109, 135} = element_bridge+1 Using Progression 1;
//+
Transfinite Curve {158, 124, 140, 163, 133, 134} = 4 Using Progression 1;
//+
Transfinite Surface {14};
//+
Transfinite Surface {17};
//+
Transfinite Surface {20};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3};
//+
Transfinite Surface {4};
//+
Transfinite Surface {27};
//+
Transfinite Surface {22};
//+
Transfinite Surface {5};
//+
Transfinite Surface {1};
//+
Transfinite Surface {26};
//+
Transfinite Surface {25};
//+
Transfinite Surface {7};
//+
Transfinite Surface {6};
//+
Transfinite Surface {12};
//+
Transfinite Surface {28};
//+
Transfinite Surface {10};
//+
Transfinite Surface {19};
//+
Transfinite Surface {9};
//+
Transfinite Surface {18};
//+
Transfinite Surface {13};
//+
Physical Surface("Volume", 172) = {1, 27, 2, 14, 17, 3, 4, 20, 22, 5, 26, 6, 12, 13, 18, 28, 10, 7, 25, 9, 19};
//+
Physical Surface("TO", 173) = {1, 5};
//+
Physical Surface("Copper", 174) = {2, 3, 4, 9, 7, 6, 18};
//+
Physical Surface("Solder", 175) = {27, 22, 26, 25};
//+
Physical Surface("Ceramic", 176) = {14, 17, 20, 12, 13, 28, 10, 19};
//+
Physical Surface("SemiconductorP+", 177) = {1};
//+
Physical Surface("SemiconductorN-", 178) = {5};
//+
Physical Curve("ElectrodeMaxX", 179) = {162};
//+
Physical Curve("ElectrodeMinX", 180) = {159};
//+
Physical Curve("Ceramic_Bottom", 181) = {140, 161, 171, 131, 134};
//+
Physical Curve("Ceramic_Top", 182) = {92, 135, 142};
//+
Physical Curve("Contact_CeramicCopper", 183) = {89, 109, 111, 119, 133, 170, 157, 124};

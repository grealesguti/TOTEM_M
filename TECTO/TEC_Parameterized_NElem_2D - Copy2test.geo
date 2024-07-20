// Gmsh project created on Mon Nov 06 18:40:42 2023

//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+
//+ Geometry definition
//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+

Copper_Height=0.1;
Ceramic_Height=0.2;
Semiconductor_Height=1.2;
Semiconductor_Width=1;
Copper_bridge=0.2;
electrode_width=0.1;
solder_height=0.05;


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
Extrude {-0.1, 0, 0} {
  Curve{91};  Recombine;
}
//+
Extrude {+0.1, 0, 0} {
  Curve{141};  Recombine;
}
Extrude {-0.1, 0, 0} {
  Curve{177};  Recombine;
}
//+
Extrude {0.1, 0, 0} {
  Curve{174};  Recombine;
}

Recursive Delete {
  Surface{30};Surface{29};
}//+
Coherence;

//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+
//+ Mesh Sizing
//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+

element_copper=3;
element_ceramic=4;
element_semi_h=25;
element_semi_w=25;
element_bridge=5;
element_solder=2;

Transfinite Curve {107, 105, 113, 114} = element_semi_h+1 Using Progression 1;
//+
Transfinite Curve {164, 165, 166, 167, 146, 145, 168, 169} = element_solder+1 Using Progression 1;
//+
Transfinite Curve {172, 147, 110, 175, 152, 162, 151, 155, 154, 159} = element_copper+1 Using Progression 1;
//+
Transfinite Curve {139, 126, 160, 129, 130, 132, 182, 178, 136, 90, 176, 185} = element_ceramic+1 Using Progression 1;
//+
Transfinite Curve {177, 173, 149, 106, 104, 156, 156, 157, 161, 161, 131, 131, 119, 119, 153, 153, 116, 116, 115, 115, 103, 174, 174, 179, 179} = element_semi_w+1 Using Progression 1;
//+
Transfinite Curve {158, 124, 140, 140, 183, 183, 184, 180, 180, 181, 163, 163, 133, 134, 134} = 4 Using Progression 1;
//+
Transfinite Curve {170, 170, 171, 150, 150, 109, 109, 135, 135} = element_bridge Using Progression 1;
//+
Transfinite Surface {32};
//+
Transfinite Surface {14};
//+
Transfinite Surface {17};
//+
Transfinite Surface {20};
//+
Transfinite Surface {31};
//+
Transfinite Surface {4};
//+
Transfinite Surface {3};
//+
Transfinite Surface {2};
//+
Transfinite Surface {27};
//+
Transfinite Surface {22};
//+
Transfinite Surface {5};
//+
Transfinite Surface {1};
//+
Transfinite Surface {25};
//+
Transfinite Surface {26};
//+
Transfinite Surface {6};
//+
Transfinite Surface {18};
//+
Transfinite Surface {13};
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
Transfinite Surface {7};


//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+
//+ Physical names // Surfaces
//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+
Physical Surface("Volume", 186) = {1, 27, 2, 14, 32, 17, 3, 20, 4, 22, 31, 5, 26, 6, 12, 13, 18, 28, 10, 7, 25, 9, 19};
//+
Physical Surface("Copper", 187) = {2, 4, 3, 6, 18, 7, 9};
//+
Physical Surface("Solder", 188) = {27, 22, 26, 25};
//+
Physical Surface("Ceramic", 189) = {14, 32, 17, 20, 31, 13, 12, 28, 10, 19};
//+
Physical Surface("SemiconductorP+", 190) = {1};
//+
Physical Surface("SemiconductorN-", 191) = {5};
//+
Physical Surface("TO", 192) = {1, 5};


//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+
//+ Physical names // Edges
//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+
//+
Physical Curve("Ceramic_Top", 193) = {184, 177, 135, 179, 180};
//+
Physical Curve("Ceramic_Bottom", 194) = {140, 161, 171, 131, 134};
//+
Physical Curve("ElectrodeMinX", 195) = {159};
//+
Physical Curve("ElectrodeMaxX", 196) = {162};
//+
Physical Curve("Contact_CeramicCopper", 197) = {173, 109, 174, 119, 133, 157, 124};
//+
Physical Curve("Symm_X", 198) = {159, 139, 132, 162, 182, 185};
//+

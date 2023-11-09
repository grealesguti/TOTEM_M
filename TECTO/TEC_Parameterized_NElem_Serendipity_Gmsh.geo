// Gmsh project created on Mon Oct 16 16:55:41 2023
SetFactory("OpenCASCADE");
//+
//+

xelem=24;
xelem_sub=xelem+1;
yelem=24;
zelem=12;
zelem_sub=zelem+1;
copper_elem=3;
copper_sub=copper_elem+1;
Height_Semiconductor=1.2;
Height_Copper=0.1;
Height_Solder=0.05;
Extra_Semiconductor=1.2;
Extra_Semiconductor_el=1;
El_TopUnion=5;


Box(1) = {0, 0, 0, 1, Height_Copper, 0.5};


Extrude {0, Height_Solder, 0} {
  Surface{4}; Layers {1}; Recombine;
}
//+
Extrude {0, Height_Semiconductor, 0} {
  Surface{11}; Layers {yelem}; Recombine;
}
//+
Extrude {0, Height_Solder, 0} {
  Surface{16}; Layers {1}; Recombine;
}
//+
Transfinite Surface {2};
//+
Transfinite Surface {6};
//+
Transfinite Surface {3};
//+
Transfinite Surface {1};
//+
Transfinite Surface {21};
//+
Transfinite Surface {5};
//+
//+
Transfinite Surface {4};
//+
//+
Extrude {+0, +Height_Copper, 0} {
  Surface{21}; Layers {copper_elem}; Recombine;
}
//+
Extrude {0.2, 0, 0} {
  Surface{23}; Layers {El_TopUnion}; Recombine;
}
//+
Extrude {1, 0, 0} {
  Surface{31}; Layers {xelem}; Recombine;
}
//+
Extrude {0, -Height_Solder, 0} {
  Surface{34}; Layers {1}; Recombine;
}
//+
Extrude {0, -Height_Semiconductor, 0} {
  Surface{41}; Layers {yelem}; Recombine;
}
//+
Extrude {0, -Height_Solder, 0} {
  Surface{46}; Layers {1}; Recombine;
}
//+
Extrude {0, -Height_Copper, 0} {
  Surface{51}; Layers {copper_elem}; Recombine;
}
//+
Extrude {0.1, 0, 0} {
  Surface{55}; Layers {1}; Recombine;
}
//+
Extrude {-0.1, 0, 0} {
  Surface{1}; Layers {1}; Recombine;
}
//+
Transfinite Curve {10, 12, 19, 27, 35, 43, 39, 57, 54, 54, 53, 53, 63, 63, 71, 71, 79, 87, 31, 31, 24, 24, 23, 23, 15, 15, 11, 11, 9, 90, 90, 82, 82, 74, 66, 66, 56, 56} = xelem_sub Using Progression 1;
//+
Transfinite Curve {107, 103, 1, 3, 20, 7, 17, 17, 28, 36, 44, 41, 51, 33, 25, 75, 83, 83, 91, 91, 5, 5, 92, 92, 100, 100, 84, 84, 99, 99, 76, 76, 68, 59, 59, 60, 60, 52, 52, 67, 67} = zelem_sub Using Progression 1;
//+
Transfinite Curve {6, 2, 4, 8} = copper_sub Using Progression 1;
//+
Physical Surface("Qin", 109) = {35, 30, 26};
//+
Physical Surface("HeatSink", 110) = {56, 3};
//+
Physical Surface("ElectrodeMaxX", 111) = {61};
//+
Physical Surface("ElectrodeMinX", 112) = {66};
//+
Physical Surface("Force", 113) = {25};
//+
//+
Transfinite Volume{1};
//+
Extrude {0, 0, -0.1} {
  Surface{12}; Surface{43}; Layers {Extra_Semiconductor_el}; Recombine;
}
//+
Transfinite Volume{15};
//+
Transfinite Volume{14};
//+
Transfinite Curve {114, 111, 119, 122} = yelem+1 Using Progression 1;
//+
//+
Extrude {0.05, 0, 0} {
  Surface{13}; Layers {1}; Recombine;
}
//+
Extrude {-0.05, 0, 0} {
  Surface{44}; Layers {1}; Recombine;
}

Extrude {0, 0, -0.1} {
  Surface{83}; Layers {Extra_Semiconductor_el}; Recombine;
}
//+
Extrude {0, 0, -0.1} {
  Surface{78}; Layers {Extra_Semiconductor_el}; Recombine;
}

//+
Transfinite Curve {138, 135, 127, 130} = 10 Using Progression 1;

Physical Volume("Copper", 114) = {1, 13, 11, 12, 7, 6, 5};
//+
Physical Volume("SemiconductorN-", 115) = {3,14,16,19};
//+
Physical Volume("SemiconductorP+", 116) = {9,15,17,18};
//+
Physical Volume("Solder", 117) = {2, 10, 8, 4};
//+
//Physical Volume("Volume", 118) = {13, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,14,15,16,17,18,19};
//+
Physical Volume("TO", 119) = {3, 9,14,15,16,17,18,19};
//+//+

//+
Coherence;
//+
Coherence;
//+
Extrude {0, 0.2, 0} {
  Surface{26}; Surface{30}; Surface{35}; Layers {5}; Recombine;
}
//+
Extrude {0, -0.2, 0} {
  Surface{3}; Surface{62}; Surface{56}; Surface{60}; Layers {5}; Recombine;
}
//+

Coherence;
//+
Coherence;
//+
Physical Volume("Ceramic", 228) = {24, 23,  25, 26, 22, 21, 20};
Physical Volume("Volume", 118) = {13, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,14,15,16,17,18,19,20,21,22,23,24,25,26};

//+
Physical Surface("Ceramic_Top", 229) = {109, 113, 117};
//+
Physical Surface("Ceramic_Bottom", 230) = {122, 131, 110};
//+
//+
Physical Surface("Top_Semiconductors", 231) = {41, 16};
//+
Physical Surface("top", 232) = {109, 113, 117};
//+
Physical Surface("bottom", 233) = {126, 122, 131, 135};
//+
Physical Surface("symmetry", 234) = {14, 42, 32, 114, 110, 27, 24, 107, 19, 77, 82, 37, 47, 52, 127, 132, 57, 6, 9, 120, 123, 63};

//+
Physical Surface("Ceramic_Contact", 235) = {26, 35, 3, 56, 62, 60};
//+
Physical Volume("Al2O3", 236) = {23, 24, 25, 26, 22, 21, 20};
//+
Physical Surface("Contact_CeramicTop", 237) = {35, 26, 30};
//+
Physical Surface("Contact_CeramicBottom", 238) = {3, 56, 62, 60};
//+
Physical Surface("Surface_2D", 239) = {14, 19, 24, 107, 114, 110, 27, 32, 37, 42, 47, 9, 6, 120, 123, 63, 52, 127, 132, 57};
//+
Physical Surface("Solder_2D", 240) = {9, 47, 37, 19};
//+
Physical Surface("Copper_2D", 241) = {24, 32, 27, 6, 63, 52, 57};
//+
Physical Surface("Ceramic_2D", 242) = {120, 123, 127, 132, 107, 110, 114};
//+
Physical Surface("SemiconductorL_2D", 243) = {14};
//+
Physical Surface("SemiconductorR_2D", 244) = {42};
//+
Physical Curve("VoltageL_2D", 245) = {105};
//+
Physical Curve("VoltageR_2D", 246) = {95};
//+
Physical Curve("TopCeramic_2D", 247) = {173, 176, 181};
//+
Physical Curve("BottomCeramic_2D", 248) = {194, 191, 200, 207};

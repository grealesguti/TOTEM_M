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
  Surface{4};  Recombine;
}
//+
Extrude {0, Height_Semiconductor, 0} {
  Surface{11};  Recombine;
}
//+
Extrude {0, Height_Solder, 0} {
  Surface{16};  Recombine;
}
//+
//+
Extrude {+0, +Height_Copper, 0} {
  Surface{21};  Recombine;
}
//+
Extrude {0.2, 0, 0} {
  Surface{23};  Recombine;
}
//+
Extrude {1, 0, 0} {
  Surface{31};  Recombine;
}
//+
Extrude {0, -Height_Solder, 0} {
  Surface{34};  Recombine;
}
//+
Extrude {0, -Height_Semiconductor, 0} {
  Surface{41};  Recombine;
}
//+
Extrude {0, -Height_Solder, 0} {
  Surface{46};  Recombine;
}
//+
Extrude {0, -Height_Copper, 0} {
  Surface{51}; Recombine;
}
//+
Extrude {0.1, 0, 0} {
  Surface{55};  Recombine;
}
//+
Extrude {-0.1, 0, 0} {
  Surface{1};  Recombine;
}
//+

//+
Extrude {0, 0, -0.1} {
  Surface{12}; Surface{43};  Recombine;
}

//+
Extrude {0.05, 0, 0} {
  Surface{13};  Recombine;
}
//+
Extrude {-0.05, 0, 0} {
  Surface{44}; Recombine;
}

Extrude {0, 0, -0.1} {
  Surface{83};  Recombine;
}
//+
Extrude {0, 0, -0.1} {
  Surface{78};  Recombine;
}


//+
Coherence;
//+
Coherence;
//+
Extrude {0, 0.2, 0} {
  Surface{26}; Surface{30}; Surface{35};  Recombine;
}
//+
Extrude {0, -0.2, 0} {
  Surface{3}; Surface{62}; Surface{56}; Surface{60};  Recombine;
}
//+

Coherence;
//+
Coherence;
//+
Extrude {0.2, 0, 0} {
  Surface{119}; 
}
//+
Recursive Delete {
  Volume{16}; 
}
//+
Recursive Delete {
  Curve{135}; Volume{17}; 
}
//+
Recursive Delete {
  Volume{19}; 
}
//+
Recursive Delete {
  Volume{18}; 
}
//+
Extrude {-0.1, 0, 0} {
  Surface{108}; Surface{25}; 
}
//+
Extrude {0.1, 0, 0} {
  Surface{116}; Surface{36}; 
}
Recursive Delete {
  Volume{15}; 
}
//+
Recursive Delete {
  Volume{14}; 
}

Coherence;



//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+
//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+
Transfinite Surface {43};
//+
Transfinite Surface {12};
//+
Transfinite Surface {42};
//+
Transfinite Surface {14};
//+
Transfinite Surface {15};
//+
Transfinite Surface {13};
//+
Transfinite Surface {44};
//+
Transfinite Surface {45};
//+
Transfinite Surface {45};
//+
Transfinite Surface {42};
//+
Transfinite Surface {37};
//+
Transfinite Surface {19};
//+
Transfinite Surface {24};
//+
Transfinite Surface {27};
//+
Transfinite Surface {32};
//+
Transfinite Surface {155};
//+
Transfinite Surface {150};
//+
Transfinite Surface {114};
//+
Transfinite Surface {110};
//+
Transfinite Surface {107};
//+
Transfinite Surface {141};
//+
Transfinite Surface {146};
//+
Transfinite Surface {167};
//+
Transfinite Surface {163};
//+
//+
Transfinite Surface {47};
//+
Transfinite Surface {9};
//+
Transfinite Surface {9};
//+
Transfinite Surface {6};
//+
Transfinite Surface {63};
//+
Transfinite Surface {123};
//+
//+
Transfinite Surface {159};
//+
Transfinite Surface {57};
//+
Transfinite Surface {132};
//+
Transfinite Surface {123};
//+
Transfinite Surface {63};
//+
Transfinite Surface {66};
//+
Transfinite Surface {125};
//+
Transfinite Surface {10};
//+
Transfinite Surface {20};
//+
Transfinite Surface {149};
//+
Transfinite Surface {145};
//+
Transfinite Surface {109};
//+
Transfinite Surface {117};
//+
Transfinite Surface {26};
//+
Transfinite Surface {21};
//+
Transfinite Surface {19};
//+
Transfinite Surface {20};
//+
Transfinite Surface {145};
//+
Transfinite Surface {149};
//+
Transfinite Surface {120};
//+
Transfinite Surface {163};
//+
Transfinite Surface {115};
//+
Transfinite Surface {105};
//+
Transfinite Surface {111};
//+
Transfinite Surface {113};
//+
Transfinite Surface {144};
//+
Transfinite Surface {142};
//+
Transfinite Surface {147};
//+
Transfinite Surface {149};
//+
Transfinite Surface {16};
//+
Transfinite Surface {38};
//+
Transfinite Surface {33};
//+
Transfinite Surface {35};
//+
Transfinite Surface {38};
//+
Transfinite Surface {37};
//+
Transfinite Surface {42};
//+
Transfinite Surface {32};
//+
Transfinite Surface {151};
//+
Transfinite Surface {156};
//+
Transfinite Surface {41};
//+
Transfinite Surface {38};
//+
Transfinite Surface {34};
//+
Transfinite Surface {160};
//+
Transfinite Surface {58};
//+
Transfinite Surface {133};
//+
Transfinite Surface {164};
//+
Transfinite Surface {163};
//+
Transfinite Surface {166};
//+
Transfinite Surface {122};
//+
Transfinite Surface {118};
//+
Transfinite Surface {5};
//+
Transfinite Surface {3};
//+
Transfinite Surface {7};
//+
Transfinite Surface {11};
//+
Transfinite Surface {46};
//+
Transfinite Surface {48};
//+
Transfinite Surface {160};
//+
Transfinite Surface {162};
//+
Transfinite Surface {159};
//+
Transfinite Surface {51};
//+
Transfinite Surface {4};
//+
Transfinite Surface {125};
//+
Transfinite Surface {66};
//+
Transfinite Surface {121};
//+
Transfinite Surface {168};
//+
Transfinite Surface {169};
//+
Transfinite Surface {165};
//+
Transfinite Surface {119};
//+
Transfinite Surface {2};
//+
Transfinite Surface {65};
//+
Transfinite Surface {1};
//+
Transfinite Surface {124};
//+
Transfinite Surface {126};
//+
Transfinite Surface {64};
//+
Transfinite Surface {170};
//+
Transfinite Surface {135};
//+
Transfinite Surface {134};
//+
Transfinite Surface {130};
//+
Transfinite Surface {61};
//+
Transfinite Surface {59};
//+
Transfinite Surface {154};
//+
Transfinite Surface {158};
//+
Transfinite Surface {110};
//+
Transfinite Surface {112};
//+
Transfinite Surface {111};
//+
Transfinite Surface {108};
//+
Transfinite Surface {17};
//+
Transfinite Surface {22};
//+
Transfinite Surface {22};
//+
Transfinite Surface {29};
//+
Transfinite Surface {30};
//+
Transfinite Surface {161};
//+
Transfinite Surface {49};
//+
Transfinite Surface {8};
//+
Transfinite Surface {2};
//+
Transfinite Surface {62};
//+
Transfinite Surface {55};
//+
Transfinite Surface {61};
//+
Transfinite Surface {60};
//+
Transfinite Surface {153};
//+
Transfinite Surface {152};
//+
Transfinite Surface {158};
//+
Transfinite Surface {157};
//+
Transfinite Surface {40};
//+
Transfinite Surface {116};
//+
Transfinite Surface {158};
//+
Transfinite Surface {36};
//+
Transfinite Surface {27};
//+
Transfinite Surface {18};
//+
Transfinite Surface {12};
//+
Transfinite Surface {23};
//+
Transfinite Surface {106};
//+
Transfinite Surface {12};
//+
Transfinite Surface {28};
//+
Transfinite Surface {33};
//+
Transfinite Surface {31};
//+
Transfinite Surface {39};
//+
Transfinite Surface {22};
//+
Transfinite Surface {25};
//+
Transfinite Surface {50};
//+
Transfinite Surface {8};
//+
Transfinite Surface {50};
//+
Transfinite Surface {8};
//+
Transfinite Surface {10};
//+
Transfinite Surface {49};
//+
Transfinite Surface {143};
//+
Transfinite Surface {148};
//+
Transfinite Volume{30};
//+
Transfinite Volume{31};
//+
Transfinite Volume{22};
//+
Transfinite Volume{7};
//+
Transfinite Volume{8};
//+
Transfinite Volume{9};
//+
Transfinite Volume{10};
//+
Transfinite Volume{11};
//+
Transfinite Volume{12};
//+
Transfinite Volume{26};
//+
Transfinite Volume{25};
//+
Transfinite Volume{27};
//+
Transfinite Volume{23};
//+
Transfinite Volume{1};
//+
Transfinite Volume{2};
//+
Transfinite Volume{3};
//+
Transfinite Volume{4};
//+
Transfinite Volume{5};
//+
Transfinite Volume{20};
//+
Transfinite Volume{6};
//+
Transfinite Volume{21};
//+
Transfinite Volume{29};
//+
Transfinite Volume{28};
//+
Transfinite Volume{13};
//+
Transfinite Volume{24};



//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+
//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+

Physical Volume("Volume", 259) = {3, 9, 10, 11, 25, 26, 12, 27, 23, 1, 2, 13, 24, 29, 28, 20, 5, 4, 6, 21, 22, 7, 8, 31, 30};
//+
Physical Volume("TO", 260) = {3, 9};
//+
Physical Volume("Ceramic", 261) = {28, 20, 21, 22, 30, 26, 25, 27, 23, 24};
//+
Physical Volume("Copper", 262) = {13, 1, 11, 12, 31, 7, 6, 5, 29};
//+
Physical Volume("Solder", 263) = {2, 10, 8, 4};
//+
Physical Volume("SemiconductorP+", 264) = {3};
//+
Physical Volume("SemiconductorN-", 265) = {9};
//+
Physical Surface("Contact_CeramicCopper", 266) = {26, 30, 35, 152, 143, 3, 62, 162, 60};
//+
Physical Surface("ElectrodeMinX", 267) = {66};
//+
Physical Surface("ElectrodeMaxX", 268) = {61};
//+
Physical Surface("Ceramic_Top", 269) = {109, 113, 117, 153, 144};
//+
Physical Surface("Ceramic_Bottom", 270) = {126, 122, 170, 166, 135};
//+
Physical Surface("Symmetry", 271) = {14, 42, 47, 159, 163, 6, 9, 120, 167, 132, 57, 63, 123, 19, 24, 107, 141, 146, 110, 27, 32, 37, 114, 150, 155};
//+
Physical Surface("Symmetry_X", 272) = {145, 149, 66, 125, 154, 158, 61, 134};
//+
//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+
//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+

xelem=24;
xelem_sub=xelem+1;
yelem=24;
zelem=12;
zelem_sub=zelem+1;
copper_elem=4;
cer_elem=5;

copper_sub=copper_elem+1;
Height_Semiconductor=1.2;
Height_Copper=0.1;
Height_Solder=0.05;
Extra_Semiconductor=1.2;
Extra_Semiconductor_el=1;
El_TopUnion=6;//+

Transfinite Curve {169, 173, 39, 43, 31, 35, 35, 23, 27, 27, 15, 15, 19, 19, 11, 12, 12, 9, 10, 10, 187, 191, 191, 183, 181, 57, 54, 56, 53, 66, 63, 74, 71, 82, 79, 248, 246, 253, 251, 251} = xelem Using Progression 1;
//+
Transfinite Curve {70, 73, 69, 72, 24, 21, 22, 26} = yelem Using Progression 1;
//+
Transfinite Curve {226, 226, 174, 171, 179, 225, 44, 239, 41, 52, 33, 28, 36, 231, 184, 25, 67, 60, 238, 238, 59, 244, 244, 68, 68, 20, 20, 107, 107, 3, 3, 103, 1, 17, 17, 7, 197, 192, 83, 75, 5, 249, 249, 189, 189, 76, 254, 254, 84, 84, 99, 92, 100, 100, 205, 210, 210, 51, 51} = zelem Using Progression 1;
//+
Transfinite Curve {34, 30, 29, 32, 64, 61, 65, 62, 78, 81, 77, 80, 16, 13, 18, 14} = 3 Using Progression 1;
//+
Transfinite Curve {108, 4, 105, 2, 8, 6, 247, 245, 89, 98, 86, 95, 55, 243, 58, 241, 47, 50, 40, 37, 42, 228, 38, 38, 230} = copper_elem Using Progression 1;
//+
Transfinite Curve {224, 168, 221, 221, 172, 167, 170, 177, 175, 175, 182, 237, 180, 234, 206, 199, 208, 208, 202, 202, 250, 252, 252, 188, 185, 185, 190, 190, 193, 193, 186, 195} = cer_elem Using Progression 1;
//+
Transfinite Curve {178, 176, 49, 46, 48, 45, 257, 255, 255, 258, 256} = El_TopUnion Using Progression 1;
//+
Transfinite Curve {223, 220, 223, 222, 219, 229, 227, 240, 242, 232, 235, 233, 236, 236, 96, 93, 97, 94, 209, 207, 194, 196, 101, 102, 104, 106} = 4 Using Progression 1;

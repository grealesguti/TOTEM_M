// Gmsh project created on Mon Oct 16 16:55:41 2023
SetFactory("OpenCASCADE");
//+
//+

xelem=2;
xelem_sub=xelem+1;
yelem=2;
zelem=1;
zelem_sub=zelem+1;
copper_elem=1;
copper_sub=copper_elem+1;
Height_Semiconductor=1.2;
Height_Copper=0.05;
Height_Solder=0.05;
Extra_Semiconductor=1.2;
Extra_Semiconductor_el=1;


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
  Surface{23}; Layers {1}; Recombine;
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
Physical Volume("Volume", 118) = {13, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,14,15,16,17,18,19};
//+
Physical Volume("TO", 119) = {3, 9,14,15,16,17,18,19};
//+//+

//+
Coherence;
//+
Coherence;

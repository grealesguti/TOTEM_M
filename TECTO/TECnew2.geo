// Gmsh project created on Fri Aug 18 12:35:20 2023
SetFactory("OpenCASCADE");
Merge "TECblock.step";
Coherence;
//+
SetFactory("OpenCASCADE");
//+
//+
//+ Electrode 1 Y
Extrude {0, -0.2, 0} { 
  Surface{16}; Layers {2}; Recombine;
}
//+ First Semiconductor X
Extrude {-1, 0, 0} {
  Surface{13}; Layers {10}; Recombine;
}
//+ Solder 1 Y
Extrude {0, 0.1, 0} { 
  Surface{24}; Layers {1}; Recombine;
}
//+ First Semiconductor Y
Extrude {0, 0.8, 0} {
  Surface{33}; Layers {10}; Recombine;
}
//+
Extrude {0, 0.1, 0} {
  Surface{38}; Layers {1}; Recombine;
}
//++Solder 3 Y 
Extrude {0, 0.2, 0} {
  Surface{43}; Layers {2}; Recombine;
}
//+
Extrude {0, -0.2, 0} {
  Surface{6}; Layers {2}; Recombine;
}
//+Second Semiconductor X
Extrude {1, 0, 0} {
  Surface{9}; Layers {10}; Recombine;
}
//+ Solder 4 Y 
Extrude {0, 0.1, 0} {
  Surface{56}; Layers {1}; Recombine;
}
//+ Second Semicoductor Y 
Extrude {0, 0.8, 0} {
  Surface{63}; Layers {10}; Recombine;
}
//+
Extrude {0, 0.1, 0} {
  Surface{68}; Layers {1}; Recombine;
}
//+
Extrude {0, 0.2, 0} {
  Surface{73}; Layers {2}; Recombine;
}
//+
Extrude {0.2, 0, 0} {
  Surface{46}; Layers {2}; Recombine;
}
//+
//+
Extrude {0, -0.2, 0} {
  Surface{24}; Layers {2}; Recombine;
}
//+
Extrude {0, -0.2, 0} {
  Surface{56}; Layers {2}; Recombine;
}
//+
Extrude {+0.1, 0, 0} {
  Surface{58}; Layers {2}; Recombine;
}
//+
Extrude {-0.1, 0, 0} {
  Surface{21}; Layers {2}; Recombine;
}
//+
Transfinite Curve {45, 42, 55, 55, 1, 1, 72, 80, 3, 6, 33, 79, 95, 143, 30, 30, 127, 128, 12, 9, 18, 21, 15, 15, 99, 119, 27, 24, 39, 36, 71, 120} = 9 Using Progression 1;
//+
Physical Volume("Volume", 185) = {2, 3, 4, 5, 6, 7, 14, 13, 12, 11, 10, 9, 8};
//+
Physical Volume("SemiconductorP+", 186) = {5};
//+
Physical Volume("SemiconductorN-", 187) = {11};
//+
Physical Volume("Copper", 188) = {7, 14, 13, 9, 8, 3, 2};
//+
Physical Volume("Solder", 189) = {4, 6, 12, 10};
//+
Physical Surface("Qin", 190) = {3};
//+
Physical Surface("Force", 191) = {2};
//+
Physical Surface("HeatSink", 192) = {26, 54};
//+
Physical Surface("ElectrodeXmin", 193) = {19};
//+
Physical Surface("ElectrodeXmax", 194) = {51};

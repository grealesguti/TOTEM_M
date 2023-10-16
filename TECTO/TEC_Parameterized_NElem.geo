// Gmsh project created on Fri Aug 18 12:35:20 2023
SetFactory("OpenCASCADE");
xelem = 3;
yelem = 3;
zelem = 3+1;
copper_elem = 5;
xelem_d=xelem+1;


Merge "TECblock.step";
//+ Copper -X
Extrude {-1, 0, 0} {
  Surface{13}; Layers {xelem}; Recombine;
}
//+ Electrode -X
Extrude {-0.1, 0, 0} {
  Surface{23}; Layers {1}; Recombine;
}
//+ Solder -X -Y
Extrude {0, +0.1, 0} {
  Surface{19}; Layers {1}; Recombine;
}
//+ Semiconductor -X
Extrude {0, +0.8, 0} {
  Surface{33}; Layers {yelem}; Recombine;
}
//+Solder -X +Y
Extrude {0, +0.1, 0} {
  Surface{38}; Layers {1}; Recombine;
}
//+ Copper -X +Y
Extrude {0, +0.2, 0} {
  Surface{43}; Layers {1}; Recombine;
}
//+ Copper Middle
Extrude {0.2, 0, 0} {
  Surface{46}; Layers {1}; Recombine;
}
//+ Copper +X +Y
Extrude {1, 0, 0} {
  Surface{53}; Layers {xelem}; Recombine;
}
//+ Solder +X +Y
Extrude {0, -0.1, 0} {
  Surface{56}; Layers {1}; Recombine;
}
//+ Semiconductor +X
Extrude {0, -0.8, 0} {
  Surface{63}; Layers {yelem}; Recombine;
}
//+ Solder +X -Y
Extrude {0, -0.1, 0} {
  Surface{68}; Layers {1}; Recombine;
}
//+ Copper +X -Y
Extrude {0, -0.2, 0} {
  Surface{73}; Layers {1}; Recombine;
}
//+ Electrode +X
Extrude {0.1, 0, 0} {
  Surface{77}; Layers {1}; Recombine;
}
//+
Coherence;
//+
Coherence;
//+
Coherence;
//+
Transfinite Curve {137, 145, 144, 136, 136, 72, 72, 151, 135, 71, 127, 195, 203, 204, 212, 196, 128, 211, 120, 183, 188, 176, 167, 175, 119, 159, 79, 168, 160, 80} = zelem Using Progression 1;
//+
Transfinite Surface {83};
//+
Transfinite Surface {79};

//+ MeshAdapt
//+ Delaunay
//+ Recombine
//+ None//+

Physical Volume("Volume", 213) = {3, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
//+
Physical Volume("TO", 214) = {5, 11};
//+
Physical Volume("SemiconductorN-", 215) = {5};
//+
Physical Volume("SemiconductorP+", 216) = {11};
//+
Physical Volume("Copper", 217) = {9, 8, 7, 2, 3, 13, 14};
//+
Physical Volume("Solder", 218) = {10, 6, 4, 12};
//+
Physical Surface("ElectrodeMinX", 219) = {79};
//+
Physical Surface("ElectrodeMaxX", 220) = {122};
//+
Physical Surface("HeatSink", 221) = {117, 81};
//+
Physical Surface("Qin", 222) = {106, 97, 93};
//+
Physical Surface("Force", 223) = {92};
//+
//Transfinite Curve {198, 123, 190, 201, 193, 193, 126, 126, 115, 177, 184, 184, 118, 118, 180, 180, 186, 186, 165, 157, 78, 78, 162, 162, 154, 154, 75, 75, 70, 70, 134, 134, 149, 149, 67, 67, 131, 131, 150} = xelem_d Using Progression 1;
//+
//Transfinite Surface {117};
//+
//Transfinite Surface {104};
//+
//Transfinite Surface {106};
//+
//Transfinite Surface {89};
//+
//Transfinite Surface {82};
//+
//Transfinite Surface {81};
//+
//Transfinite Surface {93};
//+
//+Coherence;//+
//+Transfinite Surface {93} Alternated;
//+Coherence;//+//+
//+Transfinite Curve {162, 165, 154, 157} = 4 Using Progression 1;
//+
//+Transfinite Surface {89} Alternated;
//+
//+Transfinite Surface {93} Alternated;

//+Coherence;//+//+
//+
//+
Extrude {0, -0.2, 0} {
  Surface{108}; Surface{93}; Layers {1}; Recombine;
}
//+
Coherence;
//+
Coherence;

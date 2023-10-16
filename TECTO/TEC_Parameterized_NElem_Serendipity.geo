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
Transfinite Surface {83};
//+
Transfinite Surface {114};
//+
Transfinite Volume{14};
//+
Transfinite Volume{13};
//+
Transfinite Volume{12};
//+
Transfinite Volume{11};
//+
Transfinite Volume{10};
//+
Transfinite Volume{9};
//+
Transfinite Volume{8};
//+
Transfinite Volume{7};
//+
Transfinite Volume{6};
//+
Transfinite Volume{5};
//+
Transfinite Volume{4};
//+
Transfinite Volume{2};
//+
Transfinite Volume{3};
//+
Transfinite Surface {113};
//+
Transfinite Surface {112};
//+
Transfinite Surface {117};
//+
Extrude {-1, 0, 0} {
  Surface{91}; 
}
//+
Coherence;
//+
Coherence;
//+
Extrude {1, 0, 0} {
  Surface{115}; 
}
//+
Coherence;
//+
Coherence;
//+
Coherence;
//+
Coherence;
//+
Recursive Delete {
  Surface{129}; 
}
//+
Recursive Delete {
  Curve{223}; 
}
//+
Recursive Delete {
  Curve{221}; 
}

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

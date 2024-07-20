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
  Surface{108}; 
}
//+
Extrude {0.1, 0, 0} {
  Surface{116}; 
}
Recursive Delete {
  Volume{15}; 
}
//+
Recursive Delete {
  Volume{14}; 
}

Coherence;




//+
//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+
//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+


//+
Translate {0, -0.05, 0} {
  Surface{63}; Volume{13}; Volume{24}; Volume{1}; Volume{23}; Volume{27}; Volume{25}; Volume{11}; Volume{12}; Volume{26}; 
}
//+
Translate {0, 0.05, 0} {
  Volume{5}; Volume{6}; Volume{21}; Volume{20}; Volume{28}; Volume{4}; Volume{8}; Volume{7}; Volume{22}; Volume{29}; 
}
//+
//+
Extrude {0, 0.05, 0} {
  Surface{196}; 
}
//+
Extrude {0, 0.05, 0} {
  Surface{175}; 
}
//+
Extrude {0.075, 0, 0} {
  Surface{13}; 
}
//+
Extrude {-0.075, 0, 0} {
  Surface{44}; 
}
//+
Extrude {0, -0.05, 0} {
  Surface{235}; 
}
//+

//+
Extrude {0, -0.05, 0} {
  Surface{241}; 
}
Coherence;//+
Extrude {0, 0, -0.2} {
  Surface{260}; Surface{284}; Surface{289}; Surface{266}; 
}
Coherence;//+//+

//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+
//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+

Physical Volume("Volume", 487) = {35, 10, 25, 23, 22, 21, 20, 13, 12, 11, 24, 9, 8, 7, 3, 4, 5, 6, 2, 34, 33, 32, 31, 26, 27, 28, 29, 30, 1, 39, 38, 37, 36};
//+
Physical Volume("TO", 488) = {35, 10, 9, 33, 32, 2, 3, 34, 39, 38, 37, 36};
//+
Physical Volume("SemiconductorP+", 489) = {32, 3, 34, 2, 36, 37};
//+
Physical Volume("SemiconductorN-", 490) = {35, 33, 9, 10, 38, 39};
//+
Physical Volume("Copper", 491) = {7, 6, 5};
//+
Physical Volume("Copper", 491) += {13, 12, 11, 1};
//+
Physical Volume("Ceramic", 492) = {29, 28, 22, 21, 20};
//+
Physical Volume("Ceramic", 492) += {24, 23, 27, 25, 26};
//+
Physical Volume("Solder", 493) = {8, 4};
//+
Physical Volume("Solder", 493) += {31, 30};
//+
Physical Surface("ElectrodeMinX", 494) = {167};
//+
Physical Surface("ElectrodeMaxX", 495) = {201};
//+
Physical Surface("Ceramic_Top", 496) = {225, 229, 221, 249, 253};
//+
Physical Surface("Ceramic_Bottom", 497) = {172, 181, 185, 191, 205};
//+
Physical Surface("Symmetry", 498) = {265, 262};
//+
Physical Surface("Symmetry_X", 499) = {254, 230, 167, 171, 201, 204};
//+
Physical Surface("Contact_CeramicCopper", 501) = {211, 215, 244};
//+
Physical Surface("Contact_CeramicCopper", 501) += {174, 190, 200, 163};

//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+
//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+//+

xelem=24;
xelem_sub=xelem+1;
yelem=24;
zelem=12;
zelem_sub=zelem+1;
copper_elem=4;
cer_elem=5;
sol_elem=3;
bridge_elem=5;


Transfinite Curve {19, 411, 403, 366, 350, 276, 274, 281, 327, 328, 362, 275, 273, 15, 407, 399, 348, 322, 321, 278, 467, 468} = xelem Using Progression 1;
//+
Transfinite Curve {423, 415, 385, 71, 291, 292, 300, 379, 370, 371, 74, 375, 381, 387, 418, 426, 374, 303, 295, 294, 485, 486} = xelem Using Progression 1;
//+
Transfinite Curve {408, 410, 413, 414, 439, 447, 406, 405, 416, 450, 417, 442, 463, 466, 473, 480, 481, 484} = yelem Using Progression 1;
//+
Transfinite Curve {323, 326, 333, 380, 382, 336, 320, 319} = copper_elem Using Progression 1;
//+
Transfinite Curve {252, 251, 270, 298, 299, 308, 311, 302, 301, 272, 260, 259} = copper_elem Using Progression 1;
//+
Transfinite Curve {363, 365, 368, 369, 360, 373, 372, 361} = sol_elem Using Progression 1;
//+
Transfinite Curve {430, 429, 435, 436, 434, 433, 432, 431} = sol_elem Using Progression 1;
//+
Transfinite Curve {400, 402, 421, 422, 397, 425, 424, 398, 455, 456, 457, 458, 453, 460, 459, 454} = 3 Using Progression 1;
//+
Transfinite Curve {340, 391, 384, 354, 349, 339, 342, 343, 347, 357, 386, 394} = cer_elem Using Progression 1;
//+
Transfinite Curve {262, 261, 279, 284, 290, 314, 265, 264, 277, 287, 293, 316} = cer_elem Using Progression 1;
//+
Physical Surface("All_surf", 500) = {223, 226, 233, 208, 262, 295, 263, 211, 229, 222, 206, 210, 209, 224, 225, 227, 228, 230, 231, 234, 235, 260, 264, 293, 296, 302, 304, 305, 217, 212, 246, 250, 242, 236, 297, 63, 283, 168, 257, 288, 281, 177, 265, 180, 182, 270, 275, 192, 197, 187, 202, 221, 219, 220, 218, 215, 207, 249, 216, 232, 213, 214, 294, 238, 286, 299, 253, 244, 290, 248, 247, 240, 254, 241, 251, 252, 269, 245, 243, 239, 237, 300, 298, 261, 287, 292, 267, 258, 282, 164, 167, 166, 11, 165, 163, 268, 259, 171, 175, 170, 284, 285, 174, 169, 172, 256, 255, 289, 291, 280, 279, 272, 277, 173, 176, 181, 194, 184, 46, 266, 179, 274, 178, 196, 186, 183, 185, 190, 273, 271, 278, 276, 199, 195, 193, 201, 198, 200, 191, 189, 188, 204, 203, 205, 307, 310, 316, 301, 308, 313, 312, 303, 315, 309, 306, 314, 311, 317, 318};
//+
Transfinite Curve {331, 332, 341, 334, 335, 344, 283, 282, 286, 285} = bridge_elem Using Progression 1;
//+
Transfinite Curve {268, 255, 258, 358, 359, 267, 254, 257, 404, 20, 412, 367, 329, 330, 351, 280, 269, 271, 401, 17, 409, 364, 324, 325, 345, 443, 444, 452, 451, 289, 288, 304, 427, 75, 419, 376, 337, 338, 346, 297, 296, 305, 428, 76, 420, 378, 377, 383, 388, 318, 313, 312, 395, 396} = 10 Using Progression 1;
//+
Transfinite Curve {250, 353, 352, 263, 249, 266, 355, 356, 256, 253, 390, 389, 306, 307, 315, 393, 392, 309, 310, 317} = 3 Using Progression 1;
//+
Transfinite Curve {438, 437, 441, 440, 472, 470, 445, 446, 448, 449, 476, 479} = 3 Using Progression 1;

Transfinite Surface {223};
Transfinite Surface {226};
Transfinite Surface {233};
Transfinite Surface {208};
Transfinite Surface {262};
Transfinite Surface {295};
Transfinite Surface {263};
Transfinite Surface {211};
Transfinite Surface {229};
Transfinite Surface {222};
Transfinite Surface {206};
Transfinite Surface {210};
Transfinite Surface {209};
Transfinite Surface {224};
Transfinite Surface {225};
Transfinite Surface {227};
Transfinite Surface {228};
Transfinite Surface {230};
Transfinite Surface {231};
Transfinite Surface {234};
Transfinite Surface {235};
Transfinite Surface {260};
Transfinite Surface {264};
Transfinite Surface {293};
Transfinite Surface {296};
Transfinite Surface {302};
Transfinite Surface {304};
Transfinite Surface {305};
Transfinite Surface {217};
Transfinite Surface {212};
Transfinite Surface {246};
Transfinite Surface {250};
Transfinite Surface {242};
Transfinite Surface {236};
Transfinite Surface {297};
Transfinite Surface {63};
Transfinite Surface {283};
Transfinite Surface {168};
Transfinite Surface {257};
Transfinite Surface {288};
Transfinite Surface {281};
Transfinite Surface {177};
Transfinite Surface {265};
Transfinite Surface {180};
Transfinite Surface {182};
Transfinite Surface {270};
Transfinite Surface {275};
Transfinite Surface {192};
Transfinite Surface {197};
Transfinite Surface {187};
Transfinite Surface {202};
Transfinite Surface {221};
Transfinite Surface {219};
Transfinite Surface {220};
Transfinite Surface {218};
Transfinite Surface {215};
Transfinite Surface {207};
Transfinite Surface {249};
Transfinite Surface {216};
Transfinite Surface {232};
Transfinite Surface {213};
Transfinite Surface {214};
Transfinite Surface {294};
Transfinite Surface {238};
Transfinite Surface {286};
Transfinite Surface {299};
Transfinite Surface {253};
Transfinite Surface {244};
Transfinite Surface {290};
Transfinite Surface {248};
Transfinite Surface {247};
Transfinite Surface {240};
Transfinite Surface {254};
Transfinite Surface {241};
Transfinite Surface {251};
Transfinite Surface {252};
Transfinite Surface {269};
Transfinite Surface {245};
Transfinite Surface {243};
Transfinite Surface {239};
Transfinite Surface {237};
Transfinite Surface {300};
Transfinite Surface {298};
Transfinite Surface {261};
Transfinite Surface {287};
Transfinite Surface {292};
Transfinite Surface {267};
Transfinite Surface {258};
Transfinite Surface {282};
Transfinite Surface {164};
Transfinite Surface {167};
Transfinite Surface {166};
Transfinite Surface {11};
Transfinite Surface {165};
Transfinite Surface {163};
Transfinite Surface {268};
Transfinite Surface {259};
Transfinite Surface {171};
Transfinite Surface {175};
Transfinite Surface {170};
Transfinite Surface {284};
Transfinite Surface {285};
Transfinite Surface {174};
Transfinite Surface {169};
Transfinite Surface {172};
Transfinite Surface {256};
Transfinite Surface {255};
Transfinite Surface {289};
Transfinite Surface {291};
Transfinite Surface {280};
Transfinite Surface {279};
Transfinite Surface {272};
Transfinite Surface {277};
Transfinite Surface {173};
Transfinite Surface {176};
Transfinite Surface {181};
Transfinite Surface {194};
Transfinite Surface {184};
Transfinite Surface {46};
Transfinite Surface {266};
Transfinite Surface {179};
Transfinite Surface {274};
Transfinite Surface {178};
Transfinite Surface {196};
Transfinite Surface {186};
Transfinite Surface {183};
Transfinite Surface {185};
Transfinite Surface {190};
Transfinite Surface {273};
Transfinite Surface {271};
Transfinite Surface {278};
Transfinite Surface {276};
Transfinite Surface {199};
Transfinite Surface {195};
Transfinite Surface {193};
Transfinite Surface {201};
Transfinite Surface {198};
Transfinite Surface {200};
Transfinite Surface {191};
Transfinite Surface {189};
Transfinite Surface {188};
Transfinite Surface {204};
Transfinite Surface {203};
Transfinite Surface {205};
Transfinite Surface {307};
Transfinite Surface {310};
Transfinite Surface {316};
Transfinite Surface {301};
Transfinite Surface {308};
Transfinite Surface {313};
Transfinite Surface {312};
Transfinite Surface {303};
Transfinite Surface {315};
Transfinite Surface {309};
Transfinite Surface {306};
Transfinite Surface {314};
Transfinite Surface {311};
Transfinite Surface {317};
Transfinite Surface {318};

Transfinite Volume {35};
Transfinite Volume {10};
Transfinite Volume {25};
Transfinite Volume {23};
Transfinite Volume {22};
Transfinite Volume {21};
Transfinite Volume {20};
Transfinite Volume {13};
Transfinite Volume {12};
Transfinite Volume {11};
Transfinite Volume {24};
Transfinite Volume {9};
Transfinite Volume {8};
Transfinite Volume {7};
Transfinite Volume {3};
Transfinite Volume {4};
Transfinite Volume {5};
Transfinite Volume {6};
Transfinite Volume {2};
Transfinite Volume {34};
Transfinite Volume {33};
Transfinite Volume {32};
Transfinite Volume {31};
Transfinite Volume {26};
Transfinite Volume {27};
Transfinite Volume {28};
Transfinite Volume {29};
Transfinite Volume {30};
Transfinite Volume {1};
Transfinite Volume {39};
Transfinite Volume {38};
Transfinite Volume {37};
Transfinite Volume {36};




//+
Recursive Delete {
  Volume{32}; 
}
//+
Recursive Delete {
  Volume{37}; 
}
//+
Recursive Delete {
  Volume{33}; 
}
//+
Recursive Delete {
  Volume{38}; 
}
//+
Transfinite Curve {465, 462, 474, 482, 483, 477, 461, 464} = 5 Using Progression 1;

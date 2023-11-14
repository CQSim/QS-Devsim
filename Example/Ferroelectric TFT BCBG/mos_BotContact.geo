device_width =     1.0;
gate_width =       0.8;

bulk_thickness	=		4e-2;
air_thickness	=		4e-1;
oxide_thickness =      2e-1;
diffusion_thickness =  1e-2;

refine_spacing 		=  1e-2;

x_bulk_left =    0.0;
x_bulk_right =   x_bulk_left + device_width;
x_center =       0.5 * (x_bulk_left + x_bulk_right);
x_gate_left =    x_center - 0.5 * (gate_width);
x_gate_right =   x_center + 0.5 * (gate_width);
x_device_left =  x_bulk_left - air_thickness;
x_device_right = x_bulk_right + air_thickness;

y_bulk_top =       0.0;
y_oxide_top =      y_bulk_top - oxide_thickness;
y_oxide_mid =      0.5 * (y_oxide_top + y_bulk_top);
y_bulk_bottom =    y_bulk_top + bulk_thickness;
y_bulk_mid =       0.5 * (y_bulk_top + y_bulk_bottom);
y_device_bottom =  y_bulk_bottom + air_thickness;
y_diffusion =      y_bulk_top + diffusion_thickness;
y_air_top =			y_bulk_bottom + air_thickness;

//// Bulk
Point(1) = {x_bulk_left,  y_bulk_bottom, 0, refine_spacing};
Point(4) = {x_bulk_right, y_bulk_bottom, 0, refine_spacing};

Point(3) = {x_gate_left,  y_bulk_top, 0, refine_spacing/6};
Point(6) = {x_gate_right, y_bulk_top, 0, refine_spacing/6};

Point(9) = {x_gate_left+refine_spacing/5,  y_bulk_top, 0, refine_spacing/2};
Point(10) = {x_gate_right-refine_spacing/5, y_bulk_top, 0, refine_spacing/2};


Point(7) = {x_bulk_right, y_bulk_top, 0, refine_spacing/2};
Point(8) = {x_bulk_left,  y_bulk_top, 0, refine_spacing/2};
//// Oxide
Point(11) = {x_bulk_right, y_oxide_top, 0, refine_spacing *1};
Point(12) = {x_bulk_left,  y_oxide_top, 0, refine_spacing *1};

Point(16) = {x_bulk_left,  y_air_top, 0, refine_spacing * 2};
Point(15) = {x_bulk_right, y_air_top, 0, refine_spacing * 2};

Line(11) = {1, 4};
Line(12) = {4, 15};
Line(13) = {15, 16};
Line(14) = {16, 1};

Line Loop(10) = {11,12,13,14};
Plane Surface(11) = {10};

Line(21) = {1, 8};
Line(22) = {8, 3};
Line(23) = {3, 9};
Line(24) = {9, 10};
Line(25) = {10,6};
Line(26) = {6, 7};
Line(27) = {7, 4};


Line Loop(20) = {21, 22, 23, 24, 25, 26, 27, -11};
Plane Surface(21) = {20};

Line(31) = {8, 12};
Line(32) = {12, 11};
Line(33) = {11, 7};

Line Loop(30) = {31, 32, 33, -26, -25, -24, -23, -22};
Plane Surface(31) = {30};



Field[1] = Attractor;
Field[1].NNodesByEdge = 200;
Field[1].EdgesList = {11};

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = refine_spacing/2;
Field[2].LcMax = refine_spacing*10;
Field[2].DistMin = refine_spacing/2;
Field[2].DistMax = refine_spacing*2;

Field[14] = Attractor;
Field[14].NNodesByEdge = 100;
Field[14].EdgesList = {22,23,25,26};

Field[15] = Threshold;
Field[15].IField = 14;
Field[15].LcMin = refine_spacing/4;
Field[15].LcMax = refine_spacing*10;
Field[15].DistMin = refine_spacing/3;
Field[15].DistMax = refine_spacing*5;


Field[3] = Attractor;
Field[3].NNodesByEdge = 2000;
Field[3].EdgesList = {24};

Field[5] = Threshold;
Field[5].IField = 3;
Field[5].LcMin = refine_spacing/3;
Field[5].LcMax = refine_spacing*10;
Field[5].DistMin = refine_spacing/10;
Field[5].DistMax = refine_spacing/5;

//Field[5] = MathEval;
//Field[5].F = Sprintf("F3^5 + %g", refine_spacing/5);


Field[7] = Min;
Field[7].FieldsList = {2,5,15};
Background Field = 7;
//Mesh.CharacteristicLengthExtendFromBoundary = 0;
//Field[4] = Threshold;
//Field[4].IField = 3;
//Field[4].LcMin = refine_spacing/2;
//Field[4].LcMax = refine_spacing*5;
//Field[4].DistMin = bulk_thickness;
//Field[4].DistMax = bulk_thickness*2;

// Box field to impose a step change in element sizes inside a box
//Field[6] = Box;
//Field[6].VIn = refine_spacing/2;
//Field[6].VOut = refine_spacing;
//Field[6].XMin = x_device_left; 
//Field[6].XMax = x_device_right;
//Field[6].YMax = y_bulk_bottom;
//Field[6].YMin = y_bulk_bottom - refine_spacing;

//Field[8] = Box;
//Field[8].VIn  = refine_spacing/2;
//Field[8].VOut = refine_spacing;
//Field[8].XMin = x_device_left;
//Field[8].XMax = x_device_right;
//Field[8].YMin = -refine_spacing*2;
//Field[8].YMax =  refine_spacing*2;


//Field[9] = Box;
//Field[9].VIn  = refine_spacing/10;
//Field[9].VOut = refine_spacing/5;
//Field[9].XMin = x_gate_left; 
//Field[9].XMax = x_gate_right;
//Field[9].YMin = y_gate_top;
//Field[9].YMax = y_bulk_top/2;

//Field[16] = Threshold;
//Field[16].IField = 9;
//Field[16].LcMin = refine_spacing/2;
//Field[16].LcMax = refine_spacing/5;
//Field[16].DistMin = refine_spacing/5;
//Field[16].DistMax = refine_spacing/2;

//// Use minimum of all the fields as the background field


// Don't extend the elements sizes from the boundary inside the domain
//Mesh.CharacteristicLengthExtendFromBoundary = 0;
//Mesh.Algorithm=8; /*Delaunay*/
//Mesh.RandomFactor=1e-5; /*perturbation*/
//Mesh.CharacteristicLengthFromPoints = 0;
//Mesh.CharacteristicLengthFromCurvature = 0;
//Mesh.CharacteristicLengthExtendFromBoundary = 0;

Physical Line("bulk_oxide_interface") = {24};
Physical Line("bulk_air_interface")   = {11};

Physical Line("body_contact")   = {13};
Physical Line("gate_contact")   = {32};
Physical Line("source_contact") = {22};
Physical Line("drain_contact")  = {26};
Physical Line("attached_source_contact") = {22};
Physical Line("attached_drain_contact")  = {26};




Physical Surface("oxide") = {31};
Physical Surface("bulk")  = {21};
Physical Surface("air")   = {11};


/////

////////+



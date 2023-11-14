device_width =     1.0;
gate_width =       0.8;

bulk_thickness	=	4e-2;
air_thickness	=	4e-1;
oxide_thickness =      2e-1;
diffusion_thickness =  1e-2;

refine_spacing 		=  1.2e-2;
surface_spacing 	=  5e-3;
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

Point(3) = {x_gate_left,  y_bulk_bottom, 0, refine_spacing};
Point(6) = {x_gate_right, y_bulk_bottom, 0, refine_spacing};

Point(9) = {x_gate_left+refine_spacing/3,  y_bulk_bottom, 0, refine_spacing};
Point(10) = {x_gate_right-refine_spacing/3, y_bulk_bottom, 0, refine_spacing};


Point(7) = {x_bulk_right, y_bulk_top, 0, refine_spacing/3};
Point(8) = {x_bulk_left,  y_bulk_top, 0, refine_spacing/3};
Point(101) = {x_bulk_right, y_bulk_top+surface_spacing, 0, refine_spacing/3};
Point(102) = {x_bulk_left,  y_bulk_top+surface_spacing, 0, refine_spacing/3};



//// Oxide

Point(301) = {x_bulk_right, y_oxide_top, 0, refine_spacing *1};
Point(302) = {x_bulk_left,  y_oxide_top, 0, refine_spacing *1};
Point(303) = {x_bulk_right, y_bulk_top-surface_spacing, 0, refine_spacing/3};
Point(304) = {x_bulk_left,  y_bulk_top-surface_spacing, 0, refine_spacing/3};
Point(16) = {x_bulk_left,  y_air_top, 0, refine_spacing * 2};
Point(15) = {x_bulk_right, y_air_top, 0, refine_spacing * 2};

Line(101) = {7,8};
Line(102) = {8,102};
Line(103) = {102,101};
Line(104) = {101,7};
Line(105) = {4,101};
Line(106) = {102,1};
Line(107) = {1,3};
Line(108) = {3,9};
Line(109) = {9,10};
Line(110) = {10,6};
Line(111) = {6,4};

Line(201) = {1,16};
Line(202) = {16,15};
Line(203) = {15,4};

Line(301) = {7,303};
Line(302) = {303,301};
Line(303) = {301,302};
Line(304) = {302,304};
Line(305) = {304,8};
Line(306) = {304,303};




Line Loop(10) = {107, 108, 109, 110, 111, -203, -202, -201};
Plane Surface(10) = {10};

Line Loop(12) = {106, 107, 108, 109, 110, 111, 105, -103};
//+
Plane Surface(12) = {12};
//+
Line Loop(13) = {102, 103, 104, 101};
//+
Plane Surface(13) = {13};
//+
Line Loop(14) = {101, -305, 306, -301};
//+
Plane Surface(14) = {14};
//+
Line Loop(15) = {304, 306, 302, 303};
//+
Plane Surface(15) = {15};



Field[1] = Attractor;
Field[1].NNodesByEdge = 200;
Field[1].EdgesList = {109};

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = refine_spacing/2;
Field[2].LcMax = refine_spacing*10;
Field[2].DistMin = refine_spacing/2;
Field[2].DistMax = refine_spacing*2;

Field[14] = Attractor;
Field[14].NNodesByEdge = 100;
Field[14].EdgesList = {107, 111};

Field[15] = Threshold;
Field[15].IField = 14;
Field[15].LcMin = refine_spacing/3;
Field[15].LcMax = refine_spacing*10;
Field[15].DistMin = refine_spacing/3;
Field[15].DistMax = refine_spacing*5;


Field[3] = Attractor;
Field[3].NNodesByEdge = 2000;
Field[3].EdgesList = {101};

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

Physical Line("bulk_oxide_interface") = {101};
Physical Line("bulk_air_interface")   = {109};

Physical Line("body_contact")   = {202};
Physical Line("gate_contact")   = {303};
Physical Line("source_contact") = {107};
Physical Line("drain_contact")  = {111};
Physical Line("attached_source_contact") = {107};
Physical Line("attached_drain_contact")  = {111};

Physical Surface("oxide") = {14,15};
Physical Surface("bulk")  = {12,13};
Physical Surface("air")   = {10};


/////

//////

//+

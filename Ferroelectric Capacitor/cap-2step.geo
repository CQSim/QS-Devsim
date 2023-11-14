body_thick	= 2.0e-1;
body_step	= 2.0e-1;
oxide_thick = 2.0e-1;
lower_with	= 2.0e-1;
higher_with	= 2.0e-1;
slope_with	= 1.0e-1;

body_lag	= 2e-2;



//Oxide middle
Point(201) = {0, 0, -0, body_lag};
Point(202) = {lower_with, 0, -0, body_lag/3};
Point(203) = {lower_with+slope_with/2, body_step/2, -0, body_lag/3};
Point(101) = {lower_with, body_step/2, -0, body_lag};
Point(102) = {lower_with+slope_with, body_step/2, -0, body_lag};
Point(204) = {lower_with+slope_with, body_step, -0, body_lag/3};
Point(205) = {lower_with+slope_with+higher_with, body_step, -0, body_lag};
Point(206) = {lower_with+slope_with+higher_with, -oxide_thick, -0, body_lag};
Point(207) = {0, -oxide_thick, -0, body_lag};

Ellipse(1) = {202, 101, 202, 203};
//+
Ellipse(2) = {204, 102, 204, 203};
//+
Line(3) = {201, 202};
//+
Line(4) = {204, 205};
//+
Line(5) = {205, 206};
//+
Line(6) = {206, 207};
//+
Line(7) = {207, 201};

Curve Loop(1) = {7, 3, 1, -2, 4, 5, 6};

Plane Surface(101) = {1}; 
//+
Physical Surface("oxide") = {101};
//+

Physical Line("top_contact") = {3, 1, 2, 4};
//+
Physical Line("bot_contact") = {6};

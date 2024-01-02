
a = 40*Pi*25/180;
b = 25.0;
n = 3;

Point(1) = {-a/2, 0.0, 0.0};
Point(2) = { a/2, 0.0, 0.0};
Point(3) = { a/2,   b, 0.0};
Point(4) = {-a/2,   b, 0.0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Curve Loop(1) = {1,2,3,4};

Plane Surface(1) = {1};

Transfinite Curve{1,2,3,4} = n+1;

Physical Curve("Γᵗ") = {1};
Physical Curve("Γʳ") = {2};
Physical Curve("Γᵇ") = {3};
Physical Curve("Γˡ") = {4};

Physical Surface("Ω") = {1};

Transfinite Surface{1};

Mesh.Algorithm = 8;
Mesh.MshFileVersion = 2;
Mesh 2;
//RecombineMesh;

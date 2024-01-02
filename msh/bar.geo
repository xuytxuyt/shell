
R = 1.0;
S = Pi*R/2;
n = 3;

Point(1) = {0, 0, 0};
Point(2) = {S, 0, 0};

Line(1) = {1,2};

Transfinite Curve{1} = n+1;
Physical Curve("Ω") = {1};
Physical Point("Γᵍ") = {1};
Physical Point("Γᵗ") = {2};
Physical Point("Γᶿ") = {2};
Physical Point("Γᴹ") = {2};
Physical Point("Γᵛ") = {2};

Mesh.Algorithm = 1;
Mesh.MshFileVersion = 2;
Mesh 1;

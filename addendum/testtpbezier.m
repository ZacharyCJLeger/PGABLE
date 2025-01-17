cp(1,1) = gapoint(0,0,0);
cp(1,2) = gapoint(0.2,1,0.5);
cp(1,3) = gapoint(0,2,0);
cp(1,4) = gapoint(0,3,0);

cp(2,1) = gapoint(1,0,-0.5);
cp(2,2) = gapoint(1,1.2,1);
cp(2,3) = gapoint(1,2,0);
cp(2,4) = gapoint(1,3,-1);

cp(3,1) = gapoint(2,0,0);
cp(3,2) = gapoint(2,1,0.5);
cp(3,3) = gapoint(2,2,0);
cp(3,4) = gapoint(3,4,-1);

light
lighting gouraud

drawtpbezier(cp,1)
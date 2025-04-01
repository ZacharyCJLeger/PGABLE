function [P0,P1,P2,P3] = planarRobot(ang1,ang2,ang3, len1,len2,len3)

GA.model(PGA);

P0 = gapoint(0,0,0);

T1 = gexp(-len1/2*e01);
T2 = gexp(-len2/2*e01);
T3 = gexp(-len3/2*e01);

L = e1^e2;
R1 = gexp(-ang1/2*L);
R2 = gexp(-ang2/2*L);
R3 = gexp(-ang3/2*L);

R = R1*T1;
P1 = R*P0*inverse(R);

R = R1*T1*R2*T2;
P2 = R*P0*inverse(R);

R = R1*T1*R2*T2*R3*T3;
P3 = R*P0*inverse(R);

clf
hold on
P1 = gapoint(0,0,0); P2 = gapoint(1,0,0); P3 = gapoint(0,1,0);
draw(P1,'r'); draw(P2,'g'); draw(P3,'b');
PGADrawPolyline([P1,P2,P3,P1],'k'); view([0 90])

L1 = P1.*e3; L1=L1/norm(L); R1 = gexp(2*pi*L1/12);
S12 = R1*P2*inverse(R1); draw(S12);

L2 = P2.*e3; L2=L2/norm(L2); R2 = gexp(2*pi*L2/12);
S23 = R2*P3*inverse(R2); draw(S23);

L3 = P3.*e3; L3=L3/norm(L3); R3 = gexp(2*pi*L3/12);
S31 = R3*P1*inverse(R3); draw(S31);

C1 = (P1+S12+P2)/3;
C2 = (P2+S23+P3)/3;
C3 = (P3+S31+P1)/3;

draw(C1,'m'); draw(C2,'m'); draw(C3,'m');
PGADrawPolyline([C1,C2,C3,C1],'r');

view([0 90])

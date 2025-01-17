clf
P1 = gapoint(0,0,0); P2 = gapoint(1,0,0); P3 = gapoint(0,1,0);
draw(P1,'r'); draw(P2,'g'); draw(P3,'b');
PGADrawPolyline([P1,P2,P3,P1],'k');

J12=join(P1,P2);
J12 = J12/norm(J12);

J23=join(P2,P3);
J23 = J23/norm(J23);

J31=join(P3,P1);
J31 = J31/norm(J31);


R = -J23/J12; Ri=inverse(R);
R2 = gexp(GAZ(glog(R))/6); R2i = inverse(R2);

R = -J31/J23;
R3 = gexp(GAZ(glog(R))/6); R3i = inverse(R3);

R = -J12/J31;
R1 = gexp(GAZ(glog(R))/6); R1i = inverse(R1);


L31 = R3i*J31*R3;
draw(L31,'c')

L13 = R1*J31*R1i;
draw(L13,'c');

EP1 = IntersectLines(L13,L31);
%draw(EP1,'c')


L23 = R2i*J23*R2;
draw(L23,'m')

L32 = R3*J23*R3i;
draw(L32,'m')

EP2 = IntersectLines(L23,L32);
%draw(EP2,'m')


L12 = R1i*J12*R1;
draw(L12,'y');

L21 = R2*J12*R2i;
draw(L21,'y');

EP3 = IntersectLines(L12,L21);
%draw(EP3);

PGADrawPolyline([EP1,EP2,EP3,EP1],'k')

view([0 90])

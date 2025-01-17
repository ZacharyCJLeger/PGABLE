clf
hold on
P1 = gapoint(1,0,0); P2 = gapoint(2,0,0); 
Q1 = gapoint(0,1,0); Q2 = gapoint(1,2,0); 

J1 = join(P1,Q2);
J2 = join(P2,Q1);
J3 = galine(1,1,1, 1,1,1);
draw(J1); draw(J2); draw(J3);

J12r = J1*J2
J12 = grade(J12r,2);
draw(J12,'r')

J31r= J3*J1
J31 = grade(J31r,2)
draw(J31,'g')

(J3.euclidean*I3).*(J31.euclidean*I3)
(J1.euclidean*I3).*(J31.euclidean*I3)

view([0 90])
clf
hold on
P1 = gapoint(1,0,0); P2 = gapoint(2,0,0); P3 = gapoint(4,0,0);
Q1 = gapoint(0,1,0); Q2 = gapoint(1,2,0); Q3 = gapoint(2,3,0);
PGADrawPolyline([P1,P3],'r'); PGADrawPolyline([Q1,Q3],'r');
PGADrawPolyline([P1,Q2],'k'); PGADrawPolyline([P1,Q3],'k');
PGADrawPolyline([P2,Q1],'k'); PGADrawPolyline([P2,Q3],'k');
PGADrawPolyline([P3,Q1],'k'); PGADrawPolyline([P3,Q2],'k');


H3 = IntersectLines(join(P1,Q2),join(P2,Q1));
H2 = IntersectLines(join(P1,Q3),join(P3,Q1));
H1 = IntersectLines(join(P2,Q3),join(P3,Q2));
draw(H1,'g');
draw(H2,'g');
draw(H3,'g')
PGADrawPolyline([H1,H3],'b');

view([0 90])

A = gapoint(0,2,0);
B = gapoint(-2*sin(2*2*pi/6),2*cos(2*2*pi/6),0);
C = gapoint(2*sin(2*2*pi/6),2*cos(2*2*pi/6),0);
draw(A); draw(B); draw(C);
A21B = 2*A/3 + 1*B/3;
A12B = 1*A/3 + 2*B/3;
draw(A21B); draw(A12B);

B21C = 2*B/3 + 1*C/3;
B12C = 1*B/3 + 2*C/3;
draw(B21C); draw(B12C);

C21A = 2*C/3 + 1*A/3;
C12A = 1*C/3 + 2*A/3;
draw(C21A); draw(C12A);

ABC = (A+B+C)/3;
draw(ABC);

PGADrawPolyline([A,B,C,A],'b','LineWidth',2);
PGADrawPolyline([A21B,B12C,C21A,A12B,B21C,C12A,A21B],'b','LineWidth',2)
axis off
view([0 90])

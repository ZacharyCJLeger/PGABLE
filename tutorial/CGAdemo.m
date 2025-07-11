function CGAdemo
% CGAdemo: a short demonstrations of the basics of CGA
%
%See also gable.

% CGABLE, Copyright (c) 2025, University of Waterloo
% Copying, use and development for non-commercial purposes permitted.
%          All rights for commercial use reserved; for more information
%          contact Stephen Mann (smann@uwaterloo.ca)
%
%          This software is unsupported.
try
clc
GA.model(CGA, true);
GAScene.clearitems();
disp('Welcome to the PGABLE demonstration of CGA!');disp(' ');
disp('This script demonstrates some of the basic features of CGA in PGABLE.');
disp('At any prompt, you may hit return to continue the demo or type');
disp('    a Matlab command.  Type ^C to exit the demo early.');
disp(' ');
disp('This demo is based on the motivating exmaple in the book Geometric');
disp('    Algebra for Computer Science, Dorst, Fontijne, and Mann.  ');
disp('    Morgan-Kaufmann, 2007.');
disp(' ');
disp(' ');
disp('Some basics to know about CGA:');
disp('    In addition to three Euclidean vectors e1,e2,e3, CGA has a point');
disp('    at the origin (no) and a point at infinity (ni).');
disp('    There are three main produce in CGA (and GA in general): an inner');
disp('    product (use .* in PGABLE), an outer product (use ^ in PGABLE),');
disp('    and a geometric product (use * in PGABLE).  Generally speaking,');
disp('    the outer product is used to combine objects to create new objects');
disp('    of higher dimension, which is how you''ll see the outer product used');
disp("    in this demo.  You'll also see the other products used in this demo,");
disp("    as well as some operation such as 'dual' that we don't explain in ");
disp("    this demo.  You should read the tutorial for a more in-depth discussion ");
disp('    of all three products and these other operations.');
disp(' ');
disp('Now, on to the demo.');

disp(' ');GAprompt;disp(' ');
disp(' ');disp(' ');disp(' ');

disp('Suppose that we have three points c1,c2,c3 in a 3D space, a line L, ');
disp('    and a plane Pi.  We would like to construct a circle C through the three ');
disp('    points, rotate C around the line L, and then reflect the whole scene in ');
disp('    the plane Pi, where we construct both L and Pi from a pair of points');
disp('    and vectors.');
disp('This demo shows how geometric algebra encodes this in its Conformal Model of');
disp('    Euclidean geometry.');
disp(' ');GAprompt;disp(' ');

c1 = gapoint(1,1,1);
c2 = gapoint(2,1,1);
c3 = gapoint(1,1.5,2);

draw(c1); draw(c2); draw(c3);
text(1,1,1.2,'c1');
text(2,1,1.2,'c2');
text(1,1.5,2.2,'c3');

a1 = gapoint(-0.5,1.5,0);
draw(a1,'r');
text(-0.5,1.5,0.2,'a1');
u = e3+0.3*e1;
len = norm(u);
PGABLEDraw.arrow(gapoint(-0.5,1.5,0,PGA),gapoint(-0.2,1.5,1,PGA));
u = u/len;
text(-0.2,1.7,1,'u');

p = gapoint(0,0,0);
draw(p,'g');
text(0,0.2,0,'p');
n = e3;
PGABLEDraw.arrow(gapoint(0,0,0,PGA),gapoint(0,0,1,PGA));
text(0,0.2,1,'n');

view([-37.5 30]); axis([-0.6000    2.2485   -0.1000    1.6000   -0.1000    2.1708 ]);
disp('In the figure, we have drawn:');
disp('   Three points c1,c2,c3, for constructing a circle;');
disp('   a point a1 and a direction u for constructing a line;');
disp('   and a point p and normal n for constructing a plane.');
disp('You may want to use Matlab to rotate the figure to see it from');
disp('   angles, both now and other times in this demo.');
disp(' ');GAprompt;disp(' ');

disp('First, we construct a circle C as the outer product of three points:');
disp(' ');
disp('>> C = c1^c2^c3; draw(C)');
C = c1^c2^c3; draw(C);

disp(' ');
disp('While we need to use an explicit draw command to draw objects, ');
disp('    we omit the draw calls in the rest of this example.');
disp(' ');GAprompt;disp(' ');
disp(' ');disp(' ');disp(' ');
disp('Next, we construct a line L through a1 in direction u:');
disp(' ');
disp('>> L = a1^u^ni; ');
L = zeroepsilons(a1^u^ni); draw(2*L);

view([-37.5 30]); axis([-0.6000    2.2485   -0.1000    1.6000   -0.1000    2.1708 ]);
disp(' ');
disp('We could have also constructed a line as the outer product of two points');
disp('    and the point at infinity (ni), since a line in CGA is really a ');
disp('    circle with infinite radius.');
disp(' ');GAprompt;disp(' ');

disp('Now we construct a plane Pi through p perpendicular to n.');
disp('We could also have construct our plane as the outer product of three ');
disp('    points and the point at infinity (ni), but instead, we construct it ');
disp('    from a point in the plane and the normal to the plane:');
disp(' ');
disp('>> Pi = p.*(n*ni); ');
Pi = p.*(n*ni); draw(4*Pi);

view([35 12]); axis([-2.4 2.6 -2 4.1 -1.8 3.3]);

disp(' ');GAprompt;disp(' ');disp(' ');disp(' ');

disp('We will now start transforming our objects.');
disp(' ');
disp('Transformations in CGA are done using rotors.  Often rotors are');
disp('    constructed using an object related to the transformation.');
disp('    For example, in our demo, we want to rotate around our line L;');
disp('    This rotor is constructed using our line L and the angle that we');
disp('    will be rotating around the line.');
disp('So we now construct this rotor R to rotate around our line L by an');
disp('    angle of pi/12 radians (15 degrees):');
disp(' ');
disp('>> R = gexp(pi/12*dual(L)/2)');
R = gexp(pi/6*dual(L)/2);

disp(' ');GAprompt;disp(' ');

disp('We apply a rotor to an object by multiplying the object by the rotor');
disp('    on the left and by the inverse of the rotor on the right.');
disp("Let's apply our rotor R to our circle C:");
disp(' ');
disp('>> C1 = R*C*inverse(R); ');
C1=zeroepsilons(R*C*inverse(R)); draw(C1);

view([35 12]); axis([-2.4 2.6 -2 4.1 -1.8 3.3]);
disp(' ');
disp("In the figure, you'll now see this second circle, C2.");
disp(' ');GAprompt;disp(' ');

disp('We now repeatedly apply R to the resulting circles:');
disp(' ');
disp('>> C2 = R*C1*inverse(R); ');
disp('>> C3 = R*C2*inverse(R); ');
disp('>> C4 = R*C3*inverse(R); ');
disp('>> C5 = R*C4*inverse(R); ');
disp('>> C6 = R*C5*inverse(R); ');

C2=zeroepsilons(R*C1*inverse(R)); draw(C2,'r');view([35 12]);
C3=zeroepsilons(R*C2*inverse(R)); draw(C3);view([35 12]);
C4=zeroepsilons(R*C3*inverse(R)); draw(C4,'k');view([35 12]);
C5=zeroepsilons(R*C4*inverse(R)); draw(C5);view([35 12]);
C6=zeroepsilons(R*C5*inverse(R)); draw(C6);view([35 12]);  axis([-2.4 2.6 -2 4.1 -1.8 3.3]);

disp(' ');
disp('You may wish to rotate the figure to look down the line L to better');
disp('    see the rotation of the circles around L.');
disp('We have drawn two of the circles in a different color to make it ');
disp('    easier to see the correspondence after reflection in the next step.');
disp(' ');GAprompt;disp(' ');

disp(' ');disp(' ');disp(' ');
disp('Finally, we will reflect the line L and the circles in the plane Pi.');
disp('In this case, the plane Pi itself can be used directly as a rotor ');
disp('    to reflect other objects through it.');
disp(' ');
disp('We first reflect the line...');
disp(' ');

disp('>> Lr = Pi*L*inverse(Pi)');
Lr = Pi*L*inverse(Pi); draw(2*Lr,'c')
view([11 8]);
disp(' ');GAprompt;disp(' ');
disp('...and now the circles:');
disp(' ');

disp('>> Cr = Pi*C*inverse(Pi);');
disp('>> C1r= Pi*C1*inverse(Pi); ');
disp('>> C2r= Pi*C2*inverse(Pi); ');
disp('>> C3r= Pi*C3*inverse(Pi); ');
disp('>> C4r= Pi*C4*inverse(Pi); ');
disp('>> C5r= Pi*C5*inverse(Pi); ');
disp('>> C6r= Pi*C6*inverse(Pi); ');

Cr=zeroepsilons(Pi*C*inverse(Pi)); draw(Cr);view([11 8]);
C1r=zeroepsilons(Pi*C1*inverse(Pi)); draw(C1r);view([11 8]);
C2r=zeroepsilons(Pi*C2*inverse(Pi)); draw(C2r,'r');view([11 8]);
C3r=zeroepsilons(Pi*C3*inverse(Pi)); draw(C3r);view([11 8]);
C4r=zeroepsilons(Pi*C4*inverse(Pi)); draw(C4r,'k');view([11 8]);
C5r=zeroepsilons(Pi*C5*inverse(Pi)); draw(C5r);view([11 8]);
C6r=zeroepsilons(Pi*C6*inverse(Pi)); draw(C6r);view([11 8]);

disp(' ');
disp('You should rotate the figure to better see the reflection of the objects.');
disp(' ');GAprompt;disp(' ');

disp('This is the end of our brief introduction to CGA in PGABLE, showing.');
disp('    points, lines, planes, and circles.  ');
disp('You should see the tutorial for a discussion of spheres and a deeper');
disp('    discussion of transformations and other operations in CGA in general.');
catch ; end

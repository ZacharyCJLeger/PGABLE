function PGAdemo
% PGAdemo: a short demonstrations of the basics of PGA
%
%See also gable.

% PGABLE, Copyright (c) 2024, University of Waterloo
% Copying, use and development for non-commercial purposes permitted.
%          All rights for commercial use reserved; for more information
%          contact Stephen Mann (smann@uwaterloo.ca)
%
%          This software is unsupported.
try
clc
GA.model(PGA);
disp('Welcome to the PGABLE demonstration of PGA!');disp(' ');
disp('This script demonstrates some of the basic features of PGA in PGABLE.');
disp('At any prompt, you may hit return to continue the demo or type');
disp('a Matlab command.  Type ^C to exit the demo early.');
disp(' ');
disp('This brief demonstration is intended to introduce you to the');
disp('the basic geometric elements of Plane-based Geometric Algebra (planes,');
disp('lines, and points) as well as lines at infinity and points at infinity.');
disp('For a more detailed discussion of PGA, see our tutorial.');
disp('It is assumed that you have run GAdemo before running PGAdemo; if you');
disp('haven''t done so, please run GAdemo now, and then return to PGAdemo.');
disp('(The important background that you need from GAdemo are the outer, inner,');
disp('and geometric products.)');

disp(' ');GAprompt;disp(' ');

disp('PGABLE has an implementation of PGA, the Geometric Algebra (R3,0,1) with');
disp('graphical representations of the geometrical objects.  PGA differs from');
disp('OGA in two important ways.  First, while the basis for PGA contains ');
disp('the basis elements of OGA, e1,e2,e3, in addition, in PGA there is');
disp('a fourth basis vector, e0, that is a null vector, meaning that ');
disp('e0.e0=0 even though e0 is not the zero vector.  You can test these');
disp('properties of e0 in PGA now by typing   inner(e0,e0)  and   e0==0');

disp(' ');GAprompt;disp(' ');

disp('The second important difference between OGA and PGA is in what the');
disp('vectors, bivectors, and trivectors represent.  A vector in PGA');
disp('represents a plane (thus the name "Plane-based Geometric Algebra);');
disp('a bivector in PGA represents a line, and a trivector in PGA represents');
disp('a point.');

disp(' ');

disp('In PGA, a vector that is a linear combination of e1,e2,e3 only is');
disp('know as a Euclidean vector.  If we have a unit Euclidean vector nv');
disp('then n = nv - d e0 represents the plane with normal n that is a');
disp('distance d from the origin.');

disp(' ');GAprompt;disp(' ');

disp('Let''s first draw the three coordinate planes, each offset from');
disp('the origin by different amounts:');
disp(' >> clf;');
clf;
disp(' >> n1=e1 - e0; draw(n1,''r'');');
n1 = e1-e0;
draw(n1,'r');
drawnow
disp(' >> n2=e2-0.5*e0; draw(n2,''m'');');
n2 = e2-0.5*e0;
draw(n2,'m');
drawnow
disp(' >> n3=e3-0.25*e0; draw(n3,''b'');');
n3 = e3-0.25*e0;
draw(n3,'b');
drawnow;
disp(' ');
pause(2);
disp('In the graphics window, you should see three planes.  The red one');
disp('represents the plane n1 with normal e1 offset from the origin by a distance');
disp('of 1; the magenta represents the plane n2 with normal e2 offset from the');
disp('origin by a distance 0.5; and the blue one represents the plane n3 with');
disp('with normal e3 offset from the origin by a distance of 0.25.');
disp('You may wish to rotate these planes in Matlab using the Matlab rotation');
disp('widget.');

disp(' ');GAprompt;disp(' ');

clc
disp('In PGA, we use the outer product to compute the intersection of two or');
disp('more planes.  Thus, if we take the outer product of two planes, we');
disp('get their line of intersection.');
pause(1);
disp('Intersecting n1 and n2 give us L12=n1^n2:');
pause(1);
L12 = n1^n2
draw(L12,'y');

disp(' ');GAprompt;disp(' ');

disp('Note two things.  First, in the graphics window, you see a yellow line');
disp('with ''hairs'' on it; these hairs indicate the direction of the line.');
disp('In this case, line L12 lies on plane n1 (the red plane) and on plane n2');
disp('(the magenta plane).');
pause(2);
disp('The second thing to note is that in the text window, the line L12 is');
disp('a bivector that has a Euclidean part, e1^e2, and an e0 part, e0(0.5e1-e2).');
disp('The Euclidean part represents the bivector perpendicular to the ');
disp('direction of the line.');
disp('Computing the Euclidean dual (by dividing by I3) of this bivector,');
disp('(e1^e2)/I3, gives the direction vector of the line.  Go ahead and ');
disp('compute this direction now:');

disp(' ');GAprompt;disp(' ');

disp('The e0 part of a line represent the offset vector (of the line ');
disp('from the origin).  However, computing this offset vector from the line');
disp('is a bit more complex; see the tutorial for details.');

disp(' ');GAprompt;disp(' ');

clc
disp('If we take the outer product of 3 planes, we get a point.  We can take');
disp('the outer product of the three planes we''ve constructed, n1^n2^n3, to');
disp('construct a point:');
pause(2);
disp(' >> P=n1^n2^n3; draw(P,''y'');');
P = n1^n2^n3
draw(P,'y');

disp(' ');GAprompt;disp(' ');

disp('Again, note two things.  First, in the graphics window, we have drawn the');
disp('point P as a yellow octahedron.  You''ll see that P lies at the ');
disp('intersection of n1, n2, and n3.');

pause(3);disp(' ');

disp('Second, in the text window, we see that P is a trivevector with a ');
disp('Euclidean part e1^e2^e3, and an e0 part.  The Euclidean trivector');
disp('e1^e2^e3 represents the point at the origin.  In a normalized point,');
disp('the coefficient of e1^e2^e3 is 1, as in this example.');
pause(2); disp(' ');
disp('The e0 part of the point contains the location of the point.  Basically,');
disp('each e0 term has a missing basis vector; the coefficient of the term is');
disp('the coordinate of the missing basis vector.  For example, P, the point we');
disp('have created, has e0(0.25e1^e2 + 0.5*e1^e3 + e2^e3), and thus has');
disp('coordinates 0.25 e3, 0.5 e2, and 1 e1');

disp(' ');GAprompt;disp(' ');

clc
disp('If we intersect two parallel planes, we get a line at infinity, which');
disp('we refer to as a vanishing line.  Let''s clear the screen and draw two');
disp('parallel planes:'); disp(' ');
clf
n1 = e1;
n2 = e1-e0;
disp('n1 = e1; draw(n1);');
draw(n1);
n2 = e1-e0;
draw(n2);
disp('n2 = e1-e0; draw(n2);');

disp(' ');GAprompt;disp(' ');

disp('We know intersect planes n1,n2 (using the outer product) and draw the');
disp('result:'); disp(' ');
nv = n1^n2
draw(nv);
disp('nv = n1^n2'); 
disp('draw(nv);');
disp(' ');
disp('In the text window, note that the line nv has only an e0 component,');
disp('where as the line we saw earlier also had a Euclidean component.');
disp('The e0 component of the earlier line was an offset from the origin');
disp('to the closest point on the line.  For nv, the e0 component gives');
disp('the direction to the line at infinity.');
disp(' ');
disp('In the graphics window, we have chosen to draw nv as a dashed line ');
disp('around the Matlab drawing volume.  Technically, where this line would');
disp('appear depends on where you view it from; we have chosen to draw this');
disp('line as if viewed from the center of the Matlab drawing volume.');

disp(' ');GAprompt;disp(' ');

clc
disp('If we intersect a line at infinity with another plane, we get a');
disp('point at infinity, also known as a vanishing point:');
disp(' ');
disp('pv = nv^e2');
pv = nv^e2
draw(pv);
view([-15 40])

disp(' ');GAprompt;disp(' ');

disp('Comparing pv to the P we saw earlier,');
P
disp('we see that pv only has an e0 component and doesn''t have a purely Euclidean');
disp('component.  The e0 component of pv is a bivector that is perpendicular to');
disp('the direction in which pv lies at infinity.');
disp('We draw pv as a star on the boundary of the Matlab box, with a short ');
disp('line segment to indicate this direction.  As expected, pv lies on the ');
disp('line nv (since pv is the intersection of nv and a plane).');

disp(' ');GAprompt;disp(' ');
disp('This is the end of our brief introduction to PGA in PGABLE, showing.');
disp('planes, lines, points, vanishing lines, and vanishing planes.  ');
disp('Transformations in PGA are similar to OGA, but for a further discussion');
disp('on these transformations and on PGA in general, see our tutorial.');
catch ; end

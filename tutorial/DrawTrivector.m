function t = DrawTrivector(A,B,C,c)
% DrawTrivector(A,B,C): draw the parallelepiped bounded by A,B,C.  
%
%See also gable, DrawBivector.

% PGABLE, Copyright (c) 2024, University of Waterloo
% Copying, use and development for non-commercial purposes permitted.
%          All rights for commercial use reserved; for more information
%          contact Stephen Mann (smann@uwaterloo.ca)
%
%          This software is unsupported.

if nargin == 3
    c = 'y';
end

draw(A,'b');
draw(B,'g');
draw(C,'m');

X = [0 A.getx() A.getx()+B.getx() B.getx()];
Y = [0 A.gety() A.gety()+B.gety() B.gety()];
Z = [0 A.getz() A.getz()+B.getz() B.getz()];
patch(X,Y,Z,c);

X = [A.getx() A.getx()+C.getx() A.getx()+B.getx()+C.getx() A.getx()+B.getx() ];
Y = [A.gety() A.gety()+C.gety() A.gety()+B.gety()+C.gety() A.gety()+B.gety() ];
Z = [A.getz() A.getz()+C.getz() A.getz()+B.getz()+C.getz() A.getz()+B.getz() ];
patch(X,Y,Z,c);

X = [0 A.getx() A.getx()+C.getx() C.getx()];
Y = [0 A.gety() A.gety()+C.gety() C.gety()];
Z = [0 A.getz() A.getz()+C.getz() C.getz()];
patch(X,Y,Z,c);

X = [B.getx() B.getx()+C.getx() A.getx()+B.getx()+C.getx() A.getx()+B.getx() ];
Y = [B.gety() B.gety()+C.gety() A.gety()+B.gety()+C.gety() A.gety()+B.gety() ];
Z = [B.getz() B.getz()+C.getz() A.getz()+B.getz()+C.getz() A.getz()+B.getz() ];
patch(X,Y,Z,c);

X = [0 B.getx() B.getx()+C.getx() C.getx()];
Y = [0 B.gety() B.gety()+C.gety() C.gety()];
Z = [0 B.getz() B.getz()+C.getz() C.getz()];
patch(X,Y,Z,c);

X = [C.getx() A.getx()+C.getx() A.getx()+B.getx()+C.getx() C.getx()+B.getx() ];
Y = [C.gety() A.gety()+C.gety() A.gety()+B.gety()+C.gety() C.gety()+B.gety() ];
Z = [C.getz() A.getz()+C.getz() A.getz()+B.getz()+C.getz() C.getz()+B.getz() ];
patch(X,Y,Z,c);

t = A^B^C;

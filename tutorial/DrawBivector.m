function b = DrawBivector(A,B,c)
% DrawBivector(A,B,c) : draw the bivector A^B.  
%  A and B must be vectors.  c is an optional color argument.
%

% PGABLE, Copyright (c) 2024, University of Waterloo
% Copying, use and development for non-commercial purposes permitted.
%          All rights for commercial use reserved; for more information
%          contact Stephen Mann (smann@uwaterloo.ca)
%
%          This software is unsupported.


if nargin == 2
    c = 'y';
end
A = GAZ(A);
B = GAZ(B);
if ~GAisa(A,'vector') | ~GAisa(B,'vector')
    A
    B
    error('DrawBivector: A and B must both be vectors.');
end
% DrawBivector: draw the parallelogram bounded by A,B.
X = [0 A.getx() A.getx()+B.getx() B.getx()];
Y = [0 A.gety() A.gety()+B.gety() B.gety()];
Z = [0 A.getz() A.getz()+B.getz() B.getz()];
%biarrow(B,A,'g');
patch(X,Y,Z,c);
draw(A);
draw(B); %biarrow would be a better choice
b = A^B;

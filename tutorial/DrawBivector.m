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

if GA.model() ~= "OGA"
   error("DrawBivector only works in OGA.  Use GA.model(OGA) to switch to OGA.");
   return;
end

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

draw(A,'LineWidth',1.5);
draw(B,'LineWidth',1.5);
PGABLEDraw.arrow(gapoint(A.getx(),A.gety(),A.getz(),PGA),gapoint(A.getx()+B.getx(),A.gety()+B.gety(),A.getz()+B.getz(),PGA),'LineWidth',1.5,'Color','g');
b = A^B;

% DrawBivector: draw the parallelogram bounded by A,B.
Pts = {gapoint(0,0,0,PGA), gapoint(A.getx(),A.gety(),A.getz(),PGA), gapoint(A.getx()+B.getx(),A.gety()+B.gety(),A.getz()+B.getz(),PGA), gapoint(B.getx(),B.gety(),B.getz(),PGA)};
%biarrow(B,A,'g');
PGABLEDraw.patch(Pts,'FaceColor',c);


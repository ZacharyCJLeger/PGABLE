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

if GA.model() ~= "OGA"
   error("DrawTrivector only works in OGA.  Use GA.model(OGA) to switch to OGA.");
   return;
end

draw(A,'LineWidth',1.5,'Color','b');
draw(B,'LineWidth',1.5,'Color','g');
draw(C,'LineWidth',1.5,'Color','m');

t = A^B^C;

if 0
Pts = {gapoint(0,0,0,PGA), gapoint(A.getx(),A.gety(),A.getz(),PGA), gapoint(A.getx()+B.getx(),A.gety()+B.gety(),A.getz()+B.getz(),PGA), gapoint(B.getx(),B.gety(),B.getz(),PGA)};
PGABLEDraw.patch(Pts,'FaceColor',c);

Pts = {gapoint(A.getx(), A.gety(), A.getz(),PGA),...
       gapoint(A.getx()+C.getx(), A.gety()+C.gety(), A.getz()+C.getz(), PGA),...
       gapoint(A.getx()+B.getx()+C.getx(), A.gety()+B.gety()+C.gety(), A.getz()+B.getz()+C.getz(), PGA),...
       gapoint(A.getx()+B.getx(), A.gety()+B.gety(), A.getz()+B.getz(), PGA)};
PGABLEDraw.patch(Pts,'FaceColor',c);

Pts = {gapoint(0, 0, 0, PGA), ...
       gapoint( A.getx(), A.gety(), A.getz(), PGA), ...
       gapoint( A.getx()+C.getx(), A.gety()+C.gety(), A.getz()+C.getz(), PGA), ...
       gapoint( C.getx(), C.gety(), C.getz(), PGA)};
PGABLEDraw.patch(Pts,'FaceColor',c);

Pts = {gapoint(B.getx(), B.gety(), B.getz(),PGA),...
      gapoint(B.getx()+C.getx(), B.gety()+C.gety(), B.getz()+C.getz(),PGA),...
      gapoint(A.getx()+B.getx()+C.getx(), A.gety()+B.gety()+C.gety(), A.getz()+B.getz()+C.getz(),PGA),...
      gapoint(A.getx()+B.getx(), A.gety()+B.gety(), A.getz()+B.getz(),PGA)};
PGABLEDraw.patch(Pts,'FaceColor',c);

bPts = {gapoint(0, 0, 0, PGA), ...
       gapoint(B.getx(), B.gety(), B.getz(),PGA), ...
       gapoint(B.getx()+C.getx(), B.gety()+C.gety(), B.getz()+C.getz(),PGA), ...
       gapoint(C.getx(), C.gety(), C.getz(),PGA)};
PGABLEDraw.patch(Pts,'FaceColor',c);

Pts = {gapoint(C.getx(), C.gety(), C.getz(), PGA), ...
       gapoint(A.getx()+C.getx(), A.gety()+C.gety(), A.getz()+C.getz(), PGA), ...
       gapoint(A.getx()+B.getx()+C.getx(), A.gety()+B.gety()+C.gety(), A.getz()+B.getz()+C.getz(), PGA),...
       gapoint(C.getx()+B.getx() , C.gety()+B.gety(), C.getz()+B.getz(), PGA)};
PGABLEDraw.patch(Pts,'FaceColor',c);
else
  h = [];
  h = [h PGABLEDraw.plotline({A,A+B,A+B+C,A+C,A},'r')];
  h = [h PGABLEDraw.plotline({A+C,C,B+C,B,A+B},'r')];
  h = [h PGABLEDraw.plotline({B+C,A+B+C},'r')];
  GAScene.addstillitem(GASceneStillItem(t,h));
end


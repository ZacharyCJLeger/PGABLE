function t = DrawTrivector(A,B,C,varargin)
% DrawTrivector(A,B,C): draw the parallelepiped bounded by A,B,C.  
%
%See also gable, DrawBivector.

% PGABLE, Copyright (c) 2024, University of Waterloo
% Copying, use and development for non-commercial purposes permitted.
%          All rights for commercial use reserved; for more information
%          contact Stephen Mann (smann@uwaterloo.ca)
%
%          This software is unsupported.

arguments
  A OGA;
  B OGA;
  C OGA;
end
arguments (Repeating)
  varargin
end

if length(varargin)==0
    varargin = {'y'};
end

if GA.model() ~= "OGA"
   error("DrawTrivector only works in OGA.  Use GA.model(OGA) to switch to OGA.");
   return;
end

solid = 0;
if varargin{1} == "Solid"
  varargin = varargin(2:end);
  solid = 1;
end


draw(A,'LineWidth',1.5,'Color','b');
draw(B,'LineWidth',1.5,'Color','g');
draw(C,'LineWidth',1.5,'Color','m');

t = A^B^C;

if solid
   if length(varargin)==0
     varargin = {'y'};
   end
   Pts = {gapoint(0,0,0,PGA), gapoint(A.getx(),A.gety(),A.getz(),PGA), gapoint(A.getx()+B.getx(),A.gety()+B.gety(),A.getz()+B.getz(),PGA), gapoint(B.getx(),B.gety(),B.getz(),PGA)};
   PGABLEDraw.patch(Pts,'FaceColor',varargin{:});

   Pts = {gapoint(A.getx(), A.gety(), A.getz(),PGA),...
       gapoint(A.getx()+C.getx(), A.gety()+C.gety(), A.getz()+C.getz(), PGA),...
       gapoint(A.getx()+B.getx()+C.getx(), A.gety()+B.gety()+C.gety(), A.getz()+B.getz()+C.getz(), PGA),...
       gapoint(A.getx()+B.getx(), A.gety()+B.gety(), A.getz()+B.getz(), PGA)};
   PGABLEDraw.patch(Pts,'FaceColor',varargin{:});

   Pts = {gapoint(0, 0, 0, PGA), ...
       gapoint( A.getx(), A.gety(), A.getz(), PGA), ...
       gapoint( A.getx()+C.getx(), A.gety()+C.gety(), A.getz()+C.getz(), PGA), ...
       gapoint( C.getx(), C.gety(), C.getz(), PGA)};
   PGABLEDraw.patch(Pts,'FaceColor',varargin{:});

   Pts = {gapoint(B.getx(), B.gety(), B.getz(),PGA),...
      gapoint(B.getx()+C.getx(), B.gety()+C.gety(), B.getz()+C.getz(),PGA),...
      gapoint(A.getx()+B.getx()+C.getx(), A.gety()+B.gety()+C.gety(), A.getz()+B.getz()+C.getz(),PGA),...
      gapoint(A.getx()+B.getx(), A.gety()+B.gety(), A.getz()+B.getz(),PGA)};
   PGABLEDraw.patch(Pts,'FaceColor',varargin{:});

   Pts = {gapoint(0, 0, 0, PGA), ...
       gapoint(B.getx(), B.gety(), B.getz(),PGA), ...
       gapoint(B.getx()+C.getx(), B.gety()+C.gety(), B.getz()+C.getz(),PGA), ...
       gapoint(C.getx(), C.gety(), C.getz(),PGA)};
   PGABLEDraw.patch(Pts,'FaceColor',varargin{:});

   Pts = {gapoint(C.getx(), C.gety(), C.getz(), PGA), ...
       gapoint(A.getx()+C.getx(), A.gety()+C.gety(), A.getz()+C.getz(), PGA), ...
       gapoint(A.getx()+B.getx()+C.getx(), A.gety()+B.gety()+C.gety(), A.getz()+B.getz()+C.getz(), PGA),...
       gapoint(C.getx()+B.getx() , C.gety()+B.gety(), C.getz()+B.getz(), PGA)};
   PGABLEDraw.patch(Pts,'FaceColor',varargin{:});
else
  h = [];
  h = [h PGABLEDraw.plotline({A,A+B,A+B+C,A+C,A},'r')];
  h = [h PGABLEDraw.plotline({A+C,C,B+C,B,A+B},'r')];
  h = [h PGABLEDraw.plotline({B+C,A+B+C},'r')];
  GAScene.addstillitem(GASceneStillItem(t,h));
end


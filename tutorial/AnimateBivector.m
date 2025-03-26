function b = AnimateBivector(A,B)
% AnimateBivector(A,B) : animate the bivector A^B.  
%  A and B must be vectors.  %

% PGABLE, Copyright (c) 2025, University of Waterloo
% Copying, use and development for non-commercial purposes permitted.
%          All rights for commercial use reserved; for more information
%          contact Stephen Mann (smann@uwaterloo.ca)
%
%          This software is unsupported.

if GA.model() ~= "OGA"
   error("DrawBivector only works in OGA.  Use GA.model(OGA) to switch to OGA.");
   return;
end

A = GAZ(A);
B = GAZ(B);
if ~GAisa(A,'vector') | ~GAisa(B,'vector')
    A
    B
    error('DrawBivector: A and B must both be vectors.');
end

GAScene.usefigure()
pclf
draw(A);
plot3(0,0,0,'.');
plot3(1.0*A.getx(),1.0*A.gety(),1.0*A.getz(),'.');
plot3(1.0*B.getx(),1.0*B.gety(),1.0*B.getz(),'.');
plot3(1.0*(A.getx()+B.getx()),1.0*(A.gety()+B.gety()),1.0*(A.getz()+B.getz()),'.');
draw(1/10*B);
for i=1:10
  GAScene.deletestillitem(2);
  draw(i/10*B);
  drawnow
  pause(0.05)
end

draw(A);
draw(B);

for i=1:10
  GAScene.deletestillitem(4);
  GAScene.deletestillitem(3);

  h = PGABLEDraw.arrow(gapoint(i/10*A.getx(),i/10*A.gety(),i/10*A.getz(),PGA),gapoint(i/10*A.getx()+B.getx(),i/10*A.gety()+B.gety(),i/10*A.getz()+B.getz(),PGA),'LineWidth',1.5,'Color','g');
  GAScene.addstillitem(GASceneStillItem(A, h));

  Pts = {gapoint(0,0,0,PGA), gapoint(i/10*A.getx(),i/10*A.gety(),i/10*A.getz(),PGA), gapoint(i/10*A.getx()+B.getx(),i/10*A.gety()+B.gety(),i/10*A.getz()+B.getz(),PGA), gapoint(B.getx(),B.gety(),B.getz(),PGA)};

  h = PGABLEDraw.patch(Pts,'FaceColor','y');
  GAScene.addstillitem(GASceneStillItem(A^B, h));
  drawnow
  pause(0.05)
end

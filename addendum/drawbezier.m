function drawbezier(cp,cpf)
% drawbezier(cp,cpf)--draw a Bezier curve with control points cp.  If
%  cpf is 1, then draw the control points/polygon
% Example: drawbezier([gapoint(0,0,0),gapoint(1,0,0),gapoint(1,1,1)],1)
[row,col]=size(cp);
if ~exist('cpf','var')
  cpf=0;
end
if cpf
  for i=1:col
    draw(cp(i));
  end
  PGADrawPolyline(cp,'m');
  %plot3(cpx,cpy,cpz,'m-o','LineWidth',3.5,'MarkerSize',4);
end
hold on
for i=0:100
  t=i/100;
  pt(i+1) = bezier(cp,t);
end
axis equal
PGADrawPolyline(pt,'k');

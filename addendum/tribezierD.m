function [pt,tp] = tribezierD(cp,u,v,w)
% [pt,tp] = tribezierD(cp,u,v,w) -- evaluate a Triangular Bezier surface with control
%  points cp at u,v,w for position and tangent plane.  
%  cp should be a 2D triangular array of gapoints,
%  where Pn00 is stored at cp(1,1), P0n0 is stored at cp(n+1,1), and
%  where P00n is stored at cp(1,n+1).  u,v,w should be Barycentric coordinates
%  (ie, sum to 1).

[row,col]=size(cp);
if row~=col
  error('trbezier: row ~= col for cp');
  return;
end
n = row-1;
for i=1:n+1
  if i==n+1
	  tp = join(cp(1,1), join(cp(2,1),cp(1,2)));
  end
  for j=1:n+1-i
    for k=1:n+1-i-j+1
    pt = u*cp(j,k) + v*cp(j+1,k) + w*cp(j,k+1);
%    draw(pt,'b')
      cp(j,k) = pt;
    end
  end
end
pt = cp(1,1);

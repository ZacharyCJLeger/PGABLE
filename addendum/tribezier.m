function r = tribezier(cp,u,v,w)
% pt = tribezier(cp,u,v,w) -- evaluate a Triangular Bezier surface with control
%  points cp at u,v,w.  cp should be a 2D triangular array of gapoints,
%  where Pn00 is stored at cp(1,1), P0n0 is stored at cp(n+1,1), and
%  where P00n is stored at cp(1,n+1).  u,v,w should be Barycentric coordinates
%  (ie, sum to 1).

[row,col]=size(cp);
if row~=col
  error('trbezier: row ~= col for cp');
  return;
end
n = row-1;
%fprintf("i = 1 to %d\n",n);
for i=1:n
%  fprintf(" i=%d, j 1 to %d\n",i,n-i+1);
  for j=1:n+1-i
%    fprintf(" j=%d, k = 1 to %d\n",j,n-i-j+2);
    for k=1:n-i-j+2
%    fprintf("   i=%d j=%d k=%d (%d %d %d)\n",i,j,k, n-j-k+2,j-1,k-1);
    pt = u*cp(j,k) + v*cp(j+1,k) + w*cp(j,k+1);
%    draw(pt,'b')
      cp(j,k) = pt;
    end
  end
end
r = cp(1,1);

function [cpl,cpr] = subdivideBezier(cp,t)
% subdivideBezier(cp,t) -- subdivide a Bezier curve at t \in (0,1)
[row,col]=size(cp);

cpl(1) = cp(1);
cpr(col) = cp(col);

for i=1:col-1
  for j=1:col-i
      cp(j)=(1-t)*cp(j) + t*cp(j+1);
  end
  cpl(i+1) = cp(1);
  cpr(col-i) = cp(col-i);
end

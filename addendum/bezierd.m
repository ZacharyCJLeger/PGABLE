function [p,d] = bezierd(cp,t)
% [p,d]=bezier(cp,t) -- evaluate a Bezier curve with control points cp at t
% for position p and derivative d
[row,col]=size(cp);
for i=1:col-1
  if i==col-1
    a = cp(2)*I3; b = cp(1)*I3;
    d = (col-1)*(a.noneuclidean - b.noneuclidean;)
  end
  for j=1:col-i
      cp(j)=(1-t)*cp(j) + t*cp(j+1);
  end
end
p = cp(1);
function r = bezier(cp,t)
% bezier(cp,t) -- evaluate a Bezier curve with control points cp at t
[row,col]=size(cp);
for i=1:col-1
  for j=1:col-i
      cp(j)=(1-t)*cp(j) + t*cp(j+1);
  end
end
r = cp(1);
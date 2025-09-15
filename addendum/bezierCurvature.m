function k = bezierCurvature(cp,t)
%BEZIERCURVATURE--compute the curvature of a Bezier curve
[row,col]=size(cp);
if t~=1
  for i=1:col-1
    for j=1:col-i
        cp(j)=(1-t)*cp(j) + t*cp(j+1);
    end
  end
  fd = hdual((col-1)*(cp(2)-cp(1)));
  sd = hdual((col-1)*(col-2)*((cp(3)-cp(2))-(cp(2)-cp(1))));
else
  fd = hdual((col-1)*(cp(5)-cp(4)));
  sd = hdual((col-1)*(col-2)*((cp(5)-cp(4))-(cp(4)-cp(3))));
end
k = abs(double(GAZ(fd^sd*I3/norm(fd)^3).*e3));

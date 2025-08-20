function v = isCubicPHC(cp)
% ISCUBIECPHC(cp)-Test if the Bezier curve given by control point cp is
% a cubic Pythagorean Hodograph curve.
v = true;
v1 = hdual(cp(1))-hdual(cp(2));
l1 = sqrt(v1.*v1);
v2 = hdual(cp(2))-hdual(cp(3));
l2 = sqrt(v1.*v1);
v3 = hdual(cp(3))-hdual(cp(4));
l3 = sqrt(v1.*v1);

% Geometric progression?
r1 = l1/l2;
r2 = l2/l3;
if abs(double(r1-r2))>1e-14
  v = false;
end

% Angle between two pairs of vectors
if abs(double(v1.*v2/norm(v1)/norm(v2) - v2.*v3/norm(v2)/norm(v3))) > 1e-14
  v = false;
end

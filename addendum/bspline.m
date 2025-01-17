function r = bspline(n,cp,kv,t)
% bspline(n,cp,kv,t) -- evaluate a degree n B-spline curve with control points cp and knot vector kv at t.
% Example: cp = [gapoint(0,0,0), gapoint(1,0,0),gapoint(1,2,0),gapoint(1,2,1)]; p=bspline(2,cp,[0 0 0 0.5 1 1 1], 0.5)
% Notes:
%    If the number of cp = M, then M>n and size(kv)>=M+2*n
%    kv should be non-decreasing sequence of real numbers with no knot
%       being duplicated more than n times
 

% Need a similar one for t > kv(end-n)
if kv(n) > t 
  fprintf(2,'t < kv(n) %g %g\n',t,kv(n));
  return;
end
if kv(end-n) < t
  fprintf(2,'t > kv(end-n) %g %g\n',t,kv(end-n));
  return;
end

sz = size(kv);
L = sz(2);
% Compute the first knot after t
for I = n+2:(L-n+2)
  if t < kv(I)
    break;
  end
end

if I == L-n+2
  fprintf(2,'t > kv(L-n-1), %g %g\n',t,kv(L-n-1));
  return;
end

% Extract knots, control points for the appropriate interval of the curve
%  Possibly this should go in drawbspline, and bspline should only take
%  a single interval
Pt = cp(I-n-1:I-1);
Pk = kv(I-n:I+n-1);
size(Pt);
size(Pk);

  for k=1:n
    for i=1:n-k+1
      Pt(i) = (Pk(n+i)-t)/(Pk(n+i)-Pk(i+k-1)) * Pt(i) + (t-Pk(i+k-1))/(Pk(n+i)-Pk(i+k-1)) * Pt(i+1);
    end
  end

r = Pt(1);

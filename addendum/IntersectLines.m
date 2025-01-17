function P = IntersectLines(L1,L2)

%L1=L1/norm(L1);
%L2=L2/norm(L2);

L1e = euclidean(L1);
CP1=L1/(L1e/I3);

% (P-t e0 L1e)*L2=0 => t=(P*L2)/((e0 L1e)*L2)
d = noneuclidean(grade(CP1*L2,3));
e = noneuclidean(grade((e0*L1e)*L2,3));
t = d/e;
P = CP1-t*(e0*L1e);

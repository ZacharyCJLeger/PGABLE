function r = tpbezier(cp,u,v)
% pt = tpbezier(cp,u,v) -- evaluate a Tensor Product Bezier surface with control
%  points cp at u,v.  cp should be a 2D array of gapoints.

[row,col]=size(cp);
for i=1:row
  q(i) = bezier(cp(i,:),u);
end
%[r,v1] = bezier(q,v);
r = bezier(q,v);

return

% Wanted to compute tangent plane, but Matlab not a team player
for i=1:col
  qq(i) = bezier(cp(:,i),v);
end
[r,v2] = bezier(qq,u);

%n = double(v1.*e2 * v2.*e3 - v1.*e3 * v2.*e2)*e1 + double(v1.*e3 * v2.*e1 - v1.*e1 * v2.*e3)*e2 + double(v1.*e1 * v2.*e2 - v1.*e2 * v2.*e1)*e3;

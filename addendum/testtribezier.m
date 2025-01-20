clear qcp
clear lcp
qcp(1,1) = gapoint(0,0,0); % 200
qcp(2,1) = gapoint(1,0,1.00); % 110
qcp(3,1) = gapoint(2,0,0); % 020
qcp(1,2) = gapoint(0,1,1.00); % 101
qcp(2,2) = gapoint(1,1,1.00); % 011
qcp(1,3) = gapoint(0,2,0); % 002
n=2;

lcp(1,1) = gapoint(0,0,0); % 100
lcp(2,1) = gapoint(1,0,0); % 010
lcp(1,2) = gapoint(0,1,0); % 001
n=1;

n=2;
% Draw the control points
for i=1:n+1
  for j=1:n+1-i+1
    draw(qcp(i,j))
  end
end

% Draw the control net
for i=1:n+1
  for j=1:n+1-i
    va=[qcp(i,j), qcp(i+1,j), qcp(i,j+1), qcp(i,j)];
    PGADrawPolyline(va,'b','LineWidth',2);
  end
end

vc=0; % count the number of vertices
% miToLi is used to construct the connectivity matrix
% Sample the patch
s=10;
for i=0:s
  for j=0:s-i
    k=s-i-j;
    pt(j+1,k+1) = tribezier(qcp,i/s,j/s,k/s);
    vc = vc+1;
    x(vc) = double((pt(j+1,k+1).noneuclidean*I3).*e1);
    y(vc) = double((pt(j+1,k+1).noneuclidean*I3).*e2);
    z(vc) = double((pt(j+1,k+1).noneuclidean*I3).*e3);
    miToLi(j+1,k+1)=vc;
  end
end

% Construct the connectivity info T
tc=0; % how many entries in T
for i=1:s+1
  for j=1:s+1-i
    tc = tc+1;
    T(tc,1:3) = [miToLi(i,j), miToLi(i+1,j), miToLi(i,j+1)];
    if i+j<=s
      tc = tc+1;
      T(tc,1:3) = [miToLi(i+1,j+1), miToLi(i,j+1), miToLi(i+1,j)];
    end
  end
end

light
lighting gouraud

trisurf(T,x,y,z);

if 0
% Draw the triangles formed from the samples
for i=1:s+1
  for j=1:s+1-i
    v1 = pt(i,j);
    v2 = pt(i+1,j);
    v3 = pt(i,j+1);
    x = [double((v1.noneuclidean*I3).*e1), double((v2.noneuclidean*I3).*e1), double((v3.noneuclidean*I3).*e1), double((v1.noneuclidean*I3).*e1)];
    y = [double((v1.noneuclidean*I3).*e2), double((v2.noneuclidean*I3).*e2), double((v3.noneuclidean*I3).*e2), double((v1.noneuclidean*I3).*e2)];
    z = [double((v1.noneuclidean*I3).*e3), double((v2.noneuclidean*I3).*e3), double((v3.noneuclidean*I3).*e3), double((v1.noneuclidean*I3).*e3)];
    patch(x, y, z,'g')
    if i+j<=s
      v1 = pt(i+1,j+1);
      v2 = pt(i+1,j);
      v3 = pt(i,j+1);
      x = [double((v1.noneuclidean*I3).*e1), double((v2.noneuclidean*I3).*e1), double((v3.noneuclidean*I3).*e1), double((v1.noneuclidean*I3).*e1)];
      y = [double((v1.noneuclidean*I3).*e2), double((v2.noneuclidean*I3).*e2), double((v3.noneuclidean*I3).*e2), double((v1.noneuclidean*I3).*e2)];
      z = [double((v1.noneuclidean*I3).*e3), double((v2.noneuclidean*I3).*e3), double((v3.noneuclidean*I3).*e3), double((v1.noneuclidean*I3).*e3)];
      patch(x,y,z,'g')
    end
  end
end
end

[pt,tp] = tribezierD(qcp,1/3,1/3,1/3);
draw(pt,'r');
draw(2*tp,pt,'c');

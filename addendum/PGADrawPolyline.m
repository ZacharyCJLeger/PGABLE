function PGADrawPolyline(PA,varargin)
[r,col] = size(PA);
for i=1:col
  cpx(i)=double(PA(i).noneuclidean*I3.*e1);
  cpy(i)=double(PA(i).noneuclidean*I3.*e2);
  cpz(i)=double(PA(i).noneuclidean*I3.*e3);
end
plot3(cpx,cpy,cpz,varargin{:})

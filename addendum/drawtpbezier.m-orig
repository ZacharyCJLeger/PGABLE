function drawtpbezier(cp,cpf)
if ~exist('cpf','var')
  cpf=0;
end

if cpf
  [row,col]=size(cp);
  for i=1:row
    for j=1:col
      cpx(i,j)=double(cp(i,j).noneuclidean*I3.*e1);
      cpy(i,j)=double(cp(i,j).noneuclidean*I3.*e2);
      cpz(i,j)=double(cp(i,j).noneuclidean*I3.*e3);
      draw(cp(i,j));
    end
    plot3(cpx(i,:),cpy(i,:),cpz(i,:),'k-','LineWidth',1.5);
  end
  for j=1:col
    plot3(cpx(:,j),cpy(:,j),cpz(:,j),'k-','LineWidth',1.5);
  end
end
for i=0:20
  u=i/20;
  for j=0:20
    v=j/20;
%    [pt(i+1,j+1),n(i+1,j+1,1:3)]=tpbezier(cp,u,v);
    pt(i+1,j+1)=tpbezier(cp,u,v);
    x(i+1,j+1)=double(pt(i+1,j+1).noneuclidean*I3.*e1);
    y(i+1,j+1)=double(pt(i+1,j+1).noneuclidean*I3.*e2);
    z(i+1,j+1)=double(pt(i+1,j+1).noneuclidean*I3.*e3);
  end
end
%s = surf(x,y,z,'EdgeColor','none')
surf(x,y,z)

% Failed attempt to get shading
if 0
s.VertexNormalsMode
s.BackFaceLighting = 'lit'
s.FaceLighting = 'gouraud'
s.FaceColor = 'interp'

properties(s)
s.CDataMode = 'manual'
s.CData = repmat( [0.3 0.3 0.7] , 21 , 21 )
%s.CData = [0.3 0.3 0.7] ;
shading interp
lightangle(-45,30)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.3;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;

% For testing
for i=0:20
  for j=0:20
    s.VertexNormals(i+1,j+1,1)=0;
    s.VertexNormals(i+1,j+1,2)=0;
    s.VertexNormals(i+1,j+1,3)=1;
  end
end
end

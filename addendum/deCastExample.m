clf

cp = [gapoint(0,0,0),gapoint(1,1.0,0),gapoint(2,1.0,0),gapoint(3,0,0)];
cpo =cp;

for i=1:4
    cpx(i)=double(cp(i).noneuclidean*I3.*e1);
    cpy(i)=double(cp(i).noneuclidean*I3.*e2);
    cpz(i)=double(cp(i).noneuclidean*I3.*e3);
end
plot3(cpx,cpy,cpz,'b-o','LineWidth',3.5,'MarkerSize',10);
hold on;

cp(1) = 0.5*cp(1) + 0.5*cp(2);
cp(2) = 0.5*cp(2) + 0.5*cp(3);
cp(3) = 0.5*cp(3) + 0.5*cp(4);
for i=1:3
%    draw(cp(i),'g')
    cpx(i)=double(cp(i).noneuclidean*I3.*e1);
    cpy(i)=double(cp(i).noneuclidean*I3.*e2);
    cpz(i)=double(cp(i).noneuclidean*I3.*e3);
end
cpx = cpx(1:3);
cpy = cpy(1:3);
cpz = cpz(1:3);
plot3(cpx,cpy,cpz,'g-o','LineWidth',3.5,'MarkerSize',10);

cp(1) = 0.5*cp(1) + 0.5*cp(2);
cp(2) = 0.5*cp(2) + 0.5*cp(3);
for i=1:3
%    draw(cp(i),'g')
    cpx(i)=double(cp(i).noneuclidean*I3.*e1);
    cpy(i)=double(cp(i).noneuclidean*I3.*e2);
    cpz(i)=double(cp(i).noneuclidean*I3.*e3);
end
plot3(cpx,cpy,cpz,'g-o','LineWidth',3.5,'MarkerSize',10);

cpx = cpx(1:2);
cpy = cpy(1:2);
cpz = cpz(1:2);
plot3(cpx,cpy,cpz,'r-o','LineWidth',3.5,'MarkerSize',10);

cp(1) = 0.5*cp(1) + 0.5*cp(2);
cpx(1)=double(cp(1).noneuclidean*I3.*e1);
cpy(1)=double(cp(1).noneuclidean*I3.*e2);
cpz(1)=double(cp(1).noneuclidean*I3.*e3);

drawbezier(cpo,0)

cpx = cpx(1:1);
cpy = cpy(1:1);
cpz = cpz(1:1);
plot3(cpx,cpy,cpz,'k*','LineWidth',3.5,'MarkerSize',10);


view([0,90])
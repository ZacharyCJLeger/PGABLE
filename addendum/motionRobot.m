function motionRobot

len1=1; len2=1; len3=1;
ang1 = [0       pi/12   pi*2/12  pi*3/12    pi*4/12  pi*4/12 pi*4/12   pi*4/12   pi*4/12  pi*3/12  pi*2/24 -pi*1/24 -pi*3/24 -pi*4/24];
ang2 = [pi*4/6  pi*4/6  pi*3.5/6 pi*3/6     pi*2/6   pi/12  -pi*1/12  -pi*3/12  -pi*5/12 -pi*6/12 -pi*6/12 -pi*6/12 -pi*6/12 -pi*6/12];
ang3 = [pi*7/24 pi*5/24 pi*4/24 -pi*1/12   -pi*3/12 -pi*3/12 -pi*3/12 -pi*2/12  -pi*1/24  pi*1/12  pi*2/12  pi*3/12  pi*3/12  pi*1/12];

clf
for i=1:14
  subplot(4,4,i)
  hold on
  axis([-2 3 -3 3])
  grid on
  [P0,P1,P2,P3] = planarRobot(ang1(i),ang2(i),ang3(i), len1,len2,len3);
  DrawPolyline({P0,P1,P2,P3},['b','r','g']);
  %text(1.1*P3.getx(),1.1*P3.gety(),num2str(i));
  view([0 90])

    for i=1:i-1
      [P0,P1,P2,P3] = planarRobot(ang1(i),ang2(i),ang3(i), len1,len2,len3);
      for delta=1:5
        t=delta/5;
        da1=(1-t)*ang1(i)+t*ang1(i+1);
        da2=(1-t)*ang2(i)+t*ang2(i+1);
        da3=(1-t)*ang3(i)+t*ang3(i+1);
        [P0d,P1d,P2d,P3d] = planarRobot(da1,da2,da3, len1,len2,len3);
        DrawPolyline({P3,P3d},'m');
        P3 = P3d;
      end
    end
end
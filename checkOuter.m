function checkOuter
  be1 = [no,e1,e2,e3,ni];
  be2 = [no^e1,no^e2,no^e3,no^ni, e1^e2,e1^e3,e1^ni, e2^e3,e2^ni, e3^ni,...
         no^e1^e2,no^e1^e3,no^e1^ni, no^e2^e3,no^e2^ni, no^e3^ni, e1^e2^e3,e1^e2^ni, e1^e3^ni, e2^e3^ni,...
	 no^e1^e2^e3, no^e1^e2^ni, no^e1^e3^ni, no^e2^e3^ni, e1^e2^e3^ni, no^e1^e2^e3^ni];

for i=1:5
  for j=1:26
    gp = be1(i)*be2(j);
    gpr = be2(j)*be1(i);
    if ( (j>=11&&j<=20) || j==26 )
      gpr = -1*gpr;
    end
    op = be1(i)^be2(j);
    iz =(gp+gpr)-2*op;
    fprintf("%s %s: %s\n",char(be1(i)),char(be2(j)),char(iz));
  end
end
disp(' ');

for i=1:5
  for j=1:26
    gp = be2(j)*be1(i);
    gpr = be1(i)*be2(j);
    if ( (j>=11&&j<=20) || j==26 )
      gpr = -1*gpr;
    end
    op = be2(j)^be1(i);
    iz =(gp+gpr)-2*op;
    fprintf("%s %s: %s\n",char(be1(i)),char(be2(j)),char(iz));
  end
end
function testGeoProdCase2
% Some special cases first (which will repeat later)
fprintf("%s * %s = %s =? -1 ?\n", char(no^e1^e2^e3^ni), char(no^e1^e2^e3^ni),char(no^e1^e2^e3^ni * no^e1^e2^e3^ni));
fprintf("%s * %s = %s =? -1 ?\n", char(no^e1^e2^ni), char(no^e1^e2^ni),char(no^e1^e2^ni * no^e1^e2^ni));
fprintf("%s * %s = %s =? 1 ?\n", char(no^ni),char(no^ni), char( (no^ni)*(no^ni) ));
fprintf("%s * %s = %s =? -1 ?\n", char(e1^e2),char(e1^e2), char( (e1^e2)*(e1^e2) ));
fprintf("%s * %s = %s =? -1 ?\n", char(e1^e3),char(e1^e3), char( (e1^e3)*(e1^e3) ));
fprintf("%s * %s = %s =? -1 - no^ni ?\n", char(ni),char(no), char( (ni)*(no) ));
fprintf("%s * %s = %s =? -e1 + no^e1^ni ?\n", char(ni),char(no^e1), char( (ni)*(no^e1) ));
fprintf("%s * %s = %s =? -e1 - no^e1^ni ?\n", char(no^e2),char(e1^e2^ni), char( (no^e2)*(e1^e2^ni) ));
input("hit enter: ");

  be1 = [no, ni, ...
      no^e1,no^e2,no^e3,e1^ni, e2^ni, e3^ni...
      no^e1^e2, no^e1^e3, no^e2^e3, e1^e2^ni, e1^e3^ni, e2^e3^ni,...
      no^e1^e2^e3, e1^e2^e3^ni];
  be2 = be1;
  count=0;
  for i=1:16
    for j=1:16
	    count = count+1;
      fprintf("%d of %d: %s * %s = %s\n", count,16*16,char(be1(i)),char(be2(j)), char( be1(i)*be2(j) ));
      input("hit enter: ");
    end
  end

  be1 = [no^ni, no^e1^ni, no^e2^ni, no^e3^ni, no^e1^e2^ni, no^e1^e3^ni, no^e2^e3^ni, no^e1^e2^e3^ni];
  be2 = [1, e1, e2, e3, e1^e2, e1^e3, e2^e3, e1^e1^e2];
  for i=1:8
    for j=1:8
      fprintf("%s * %s = %s\n", char(be1(i)),char(be2(j)), char( be1(i)*be2(j) ));
      input("hit enter: ");
      fprintf("%s * %s = %s\n", char(be2(j)),char(be1(i)), char( be2(j)*be1(i) ));
      input("hit enter: ");
    end
  end
  
end

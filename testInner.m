function testInner
  be1 = [e1.*e1, no,e1,e2,e3,ni, no^e1,no^e2,no^e3,no^ni, e1^e2,e1^e3,e1^ni, e2^e3,e2^ni, e3^ni,...
         no^e1^e2,no^e1^e3,no^e1^ni, no^e2^e3,no^e2^ni, no^e3^ni, e1^e2^e3,e1^e2^ni, e1^e3^ni, e2^e3^ni,...
	 no^e1^e2^e3, no^e1^e2^ni, no^e1^e3^ni, no^e2^e3^ni, e1^e2^e3^ni, no^e1^e2^e3^ni];
  be2 = be1;
  % 6, 16, 32
  fprintf("%s . %s = %s = ???\n",char(no^e1), char(no^ni), char((no^e1) .* (no^ni)));
  fprintf("%s . %s = %s =? no ?\n",char(no), char(no^ni), char(no.*(no^ni)));
  fprintf("%s . %s = %s =? -no^e1 ?\n",char(no), char(no^e1^ni), char(no.*(no^e1^ni)));
  fprintf("%s . %s = %s =? e2 ?\n",char(e1), char(e1^e2), char(e1.*(e1^e2)));
  fprintf("%s . %s = %s =? -no ?\n",char(e1), char(no^e1), char(e1.*(no^e1)));
  fprintf("%s . %s = %s =? no^e1^e2 ?\n",char(no), char(no^e1^e2^ni), char(no.*(no^e1^e2^ni)));
  fprintf("%s . %s = %s =? -no^e3^e1 ?\n",char(no), char(no^e1^e3^ni), char(no.*(no^e1^e3^ni)));

  for i=1:32
    for j=1:32
      fprintf("%s . %s = %s\n", char(be1(i)),char(be2(j)), char( be1(i).*be2(j) ));
      input("hit enter: ");
    end
  end
end

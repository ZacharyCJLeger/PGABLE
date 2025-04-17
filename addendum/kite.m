function kp = kite()
% KP = KITE
%  Construct the vertices of the Kite tile
	P0 = gapoint(0,0,0);
	P1 = gapoint(sqrt(3),0,0);
	P2 = gapoint(sqrt(3),1,0);

	R60 = gexp(-pi/3*(e1^e2)/2);
	R60i = inverse(R60);
	P3 = zeroepsilons(R60*P1*R60i);

	kp = {P0,P1,P2,P3,P0};

function hv = hat()
% hv = hat()
%  Compute the outer vertices of the Hat tile [cite]
%  The routine uses the kite tile to generate the hat tile, drawing the.
%  kites as it generates them.

	kp0 = kite();
	R60 = gexp(-pi/3*(e1^e2)/2);
	R60i = gexp(pi/3*(e1^e2)/2);
	kp1 = versorbatch(R60i,kp0);
	kp2 = versorbatch(R60i,kp1);
	kp3 = versorbatch(R60i,kp2);

	DrawPolyline(kp0);
	DrawPolyline(kp1);
	DrawPolyline(kp2);
	DrawPolyline(kp3);

	l023 = join(kp0{2},kp0{3});
	kp4 = versorbatch(l023,kp0);
	kp5 = versorbatch(l023,kp1);
	DrawPolyline(kp4);
	DrawPolyline(kp5);

	l13e3 = kp1{3}.*e3;
	R120 = gexp(-2*pi/3*(l13e3)/2);
	kp6 = versorbatch(R120,kp1);
	DrawPolyline(kp6);

	l61e3 = kp6{1}.*e3;
	R120 = gexp(pi/3*(l61e3)/2);
	kp7 = versorbatch(R120,kp6);
	DrawPolyline(kp7);

	hv = {kp0{3},kp0{4},kp0{5},kp3{2},kp3{3},kp3{4},kp2{3},kp2{4},...
		kp6{1},kp7{2},kp7{3},kp7{4},kp4{1},kp4{4},kp4{3}};

function r = geoall(A,B)
	r.obj1 = OGA(A);
	r.obj2 = OGA(B);
	r.comp = inner(A,B);
	r.proj = r.comp/B;
	r.rej = A - r.proj;
	r.meet = meet(A,B);
	r.join = join(A,B);

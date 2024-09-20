function r = geoall(A,B)
	r.obj1 = A;
	r.obj2 = B;
	r.comp = inner(A,B);
	r.proj = r.comp/B;
	r.rej = A - r.proj;
	r.meet = meet(A,B);
	r.join = join(A,B);

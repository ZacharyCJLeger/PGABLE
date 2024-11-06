function L = line(l1, l2, l3, m1, m2, m3)
	%LINE -  Construct a line.
    %   Can be called by passing in 2 3D vectors or 6 scalars, as such:
	%   L = line(l1,l2,l3, m1,m2,m3);
	%   L = plucker(l,m);
    %      where l = [l1, l2, l3] and m = [m1, m2, m3].
    %   The vector [l1, l2, l3] is the direction of the line.
    %   The vector [m1, m2, m3] is the moment of the line.

	if nargin==2
		m1 = l2(1);
        m2 = l2(2);
        m3 = l2(3);

        l3 = l1(3);
		l2 = l1(2);
		l1 = l1(1);
	end
	L = point(m1, m2, m3).*(l1*e1 + l2*e2 + l3*e3);
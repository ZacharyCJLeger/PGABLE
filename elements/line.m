function L = line(l1, l2, l3, p1, p2, p3)
	%LINE -  Construct a line.
    %   Can be called by passing in 2 3D vectors or 6 scalars, as such:
	%   L = line(l1,l2,l3, p1,p2,p3);
	%   L = line(l,p);
    %      where l = [l1, l2, l3] and p = [p1, p2, p3].
    %   The vector [l1, l2, l3] is the direction of the line.
    %   The vector [p1, p2, p3] is a point on the line

	if nargin==2
		p1 = l2(1);
        p2 = l2(2);
        p3 = l2(3);

        l3 = l1(3);
		l2 = l1(2);
		l1 = l1(1);
	end
	L = point(p1, p2, p3).*(l1*e1 + l2*e2 + l3*e3);

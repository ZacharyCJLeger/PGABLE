function DrawPolyline(P,c)
% DrawPolyline(P,c) - draw a polyline along the tips of the vector arguments
%  P: a cell array of GA vectors
%  c: an optional color argument; c may either be a single color character, 
%     or c can be an array of color characters of size one less than the
%     size of P.
%

% GABLE, Copyright (c) 1999, University of Amsterdam
% Copying, use and development for non-commercial purposes permitted.
%          All rights for commercial use reserved; for more information
%          contact Leo Dorst (leo@wins.uva.nl).
%
%          This software is unsupported.

l = length(P);
if l < 2
	error('DrawPolyline: requires at least two vector arguments');
end

if nargin == 1
     for i=1:l
	     c(i) = 'b';
     end
elseif nargin == 2
     if length(c) ~= l-1
        for i=2:l
	     c(i) = c(1);
        end
     end
elseif nargin > 2
     error('DrawPolyline: takes only 2 arguments');
end
x = zeros(1,2);
y = zeros(1,2);
z = zeros(1,2);
hold on;
for i = 1:l
	v = P{i};
        if (GA.model()=="OGA" && ~isgrade(GAZ(v),1)) 
                error('DrawPolyline: all OGA objects in P must be vectors.');
	elseif (GA.model()=="PGA" && ~isgrade(GAZ(v),3))
                error('DrawPolyline: all PGA objects in P must be points.');
        end
	if i==1
		x(1) = getx(v);
		y(1) = gety(v);
		z(1) = getz(v);
	else
		x(2) = getx(v);
		y(2) = gety(v);
		z(2) = getz(v);
		plot3(x,y,z,'Color',c(i-1),"LineWidth",2);
		x(1) = x(2);
		y(1) = y(2);
		z(1) = z(2);
	end
end
axis('equal');

function createPHC(drawlines)
%CREATEPHC
% Interactively create a piecewise C1 cubic Pythagorean Hodograph curve.
%  If you give an option argument of 'true', then draw the lines used
%  in constructing the rotation, etc.

   if nargin~=1
     drawlines=false;
   end
	clf;
	GAScene.view(2)
	axis([0 10 0 10]); axis square

	 % Input 3 points
	 [x,y] = ginput(1);
	 cp(1,1) = gapoint(x,y,0);
	 draw(cp(1,1)); axis([0 10 0 10]); hold on;

	 [x,y] = ginput(1);
	 cp(1,2) = gapoint(x,y,0);
	 plot(x,y,'o');
	 draw(cp(1,2)); 
	 PGADrawPolyline(cp(1,:));
	 
 	 [x,y] = ginput(1);
	 cp(1,3) = gapoint(x,y,0);
	 plot(x,y,'o');
	 draw(cp(1,3)); 
	 PGADrawPolyline(cp(1,2:3));
%print -depsc -vector cphc3pts.eps

	 % Construct rotor for direction
	 ln1 = join(cp(1,1),cp(1,2));
	 l1 = norm(ln1);
	 ln2 = join(cp(1,2),cp(1,3));
	 l2 = norm(ln2);
 	 if drawlines
		GA.epsilon_tolerance(1e-14); % There were issues with eps=1e-15
	 	 draw(ln1,cp(1,1))
	 	 draw(ln2,cp(1,2))
%print -depsc -vector cphc2lines.eps
	 end
	 axis([0 10 0 10]); view([0 90]); hold on;

	 R = ln2/ln1/l1/l2;
	 ln3 = R*ln1*inverse(R);
	 if drawlines
	 	 draw(ln3,cp(1,3),'g')
%print -depsc -vector cphcRotLine.eps
	 end
	 T = gexp(-e0*hdual(cp(1,2)-cp(1,3))/2);
	 ln4 = T*ln3*inverse(T);
	 ln4 = ln4/norm(ln4);
	 if drawlines
	    draw(ln4,cp(1,3),'m')
%print -depsc -vector cphcRotTransLine.eps
	 end
	 cp(1,4) = cp(1,3)-(l2^2/l1)*(e0*ln4);
	 draw(cp(1,4));
	 PGADrawPolyline(cp(1,3:4));
	 ln4a = join(cp(1,3),cp(1,4));
	 l3 = norm(ln4a);
	 fprintf(1,"isCubicPHC? %d\n",isCubicPHC(cp(1,:)));
	 drawbezier(cp(1,:),1)
%print -depsc -vector cphcCurve.eps

	 nc=1;
	 domore = true;
	 while domore
	         nc = nc+1;
	 	 cp(nc,1) = cp(nc-1,4);
	 	 cp(nc,2) = 2*cp(nc-1,4)-cp(nc-1,3);
		 draw(cp(nc,2));
		 view([0 90]); hold on;
	 	 [x,y] = ginput(1);
		 cp(nc,3) = gapoint(x,y,0);
		 draw(cp(nc,3));
	 	 view([0 90]); hold on;
		 ln1 = join(cp(nc,1),cp(nc,2));
		 l1 = norm(ln1);
		 ln2 = join(cp(nc,2),cp(nc,3));
		 l2 = norm(ln2);
	
		 R = ln2/ln1/l1/l2;
		 ln3 = R*ln1*inverse(R);
		 T = gexp(-e0*hdual(cp(nc,2)-cp(nc,3))/2);
		 ln4 = T*ln3*inverse(T);
		 ln4 = ln4/norm(ln4);
		 cp(nc,4) = cp(nc,3)-(l2^2/l1)*(e0*ln4);
		 draw(cp(nc,4));
		 fprintf(1,"isCubicPHC? %d\n",isCubicPHC(cp(1,:)));
		 drawbezier(cp(nc,:),1)
		 view([0 90]); hold on;
		 reply = input("Another segment? Y/n: ",'s');
		 if isempty(reply) || reply == 'y'
		   domore=true;
		 else
		   domore=false;
		 end
	end
end

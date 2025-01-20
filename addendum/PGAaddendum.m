function PGAaddendum(GAn)
%PGAagendum: run sample code in tutorial.

% PGABLE, Copyright (c) 2025, University of Waterloo
% Copying, use and development for non-commercial purposes permitted.
%          All rights for commercial use reserved; for more information
%          contact Stephen Mann (smann@uwaterloo.ca)
%
%          This software is unsupported.

GA.model(PGA);
try
if ( GAn == 1 ) 
  GAps = 'PGAblock >> ';
     disp(">> % NAPOLEON'S THEORM");
     disp('>> clf; hold on');
     clf; hold on
     disp('>> P1 = gapoint(0,0,0); P2 = gapoint(1,0,0); P3 = gapoint(0,1,0);');
     P1 = gapoint(0,0,0); P2 = gapoint(1,0,0); P3 = gapoint(0,1,0);
     disp(">> draw(P1,'r'); draw(P2,'g'); draw(P3,'b');");
     draw(P1,'r'); draw(P2,'g'); draw(P3,'b');
     disp(">> PGADrawPolyline([P1,P2,P3,P1],'k','LineWidth',2); view([0 90])");
     PGADrawPolyline([P1,P2,P3,P1],'k','LineWidth',2); view([0 90])
     GAprompt;
     disp('>> L1 = P1.*e3; L1=L1/norm(L1); R1 = gexp(2*pi*L1/12);');
     L1 = P1.*e3; L1=L1/norm(L1); R1 = gexp(2*pi*L1/12);
     disp('>> S12 = R1*P2*inverse(R1); draw(S12);');
     S12 = R1*P2*inverse(R1); draw(S12);
     GAprompt;
     disp('>> L2 = P2.*e3; L2=L2/norm(L2); R2 = gexp(2*pi*L2/12);');
     L2 = P2.*e3; L2=L2/norm(L2); R2 = gexp(2*pi*L2/12);
     disp('>> S23 = R2*P3*inverse(R2); draw(S23);');
     S23 = R2*P3*inverse(R2); draw(S23);
     disp('>> L3 = P3.*e3; L3=L3/norm(L3); R3 = gexp(2*pi*L3/12);');
     L3 = P3.*e3; L3=L3/norm(L3); R3 = gexp(2*pi*L3/12);
     disp('>> S31 = R3*P1*inverse(R3); draw(S31);');
     S31 = R3*P1*inverse(R3); draw(S31);
     GAprompt;
     disp('>> C1 = (P1+S12+P2)/3;');
     C1 = (P1+S12+P2)/3;
     disp('>> C2 = (P2+S23+P3)/3;');
     C2 = (P2+S23+P3)/3;
     disp('>> C3 = (P3+S31+P1)/3;');
     C3 = (P3+S31+P1)/3;
     disp(">> draw(C1,'m'); draw(C2,'m'); draw(C3,'m');");
     draw(C1,'m'); draw(C2,'m'); draw(C3,'m');
     disp(">> PGADrawPolyline([C1,C2,C3,C1],'r','LineWidth',2);");
     PGADrawPolyline([C1,C2,C3,C1],'r','LineWidth',2);
     disp('>> view([0,90])');
     view([0,90])
     GAprompt;
     disp(' ');    disp('End of PGAblock sequence.  Returning to Matlab.');
elseif ( GAn == 2 ) 
  GAps = 'PGAblock >> ';
     disp(">> %PAPPUS'S THEOREM ");
     %PAPPUS'S THEOREM 
     disp('>> clf; hold on');
     clf; hold on
     disp('>> P1 = gapoint(1,0,0); P2 = gapoint(2,0,0); P3 = gapoint(4,0,0);');
     P1 = gapoint(1,0,0); P2 = gapoint(2,0,0); P3 = gapoint(4,0,0);
     disp('>> Q1 = gapoint(0,1,0); Q2 = gapoint(1,2,0); Q3 = gapoint(2,3,0);');
     Q1 = gapoint(0,1,0); Q2 = gapoint(1,2,0); Q3 = gapoint(2,3,0);
     disp(">> PGADrawPolyline([P1,P3],'r','LineWidth',2); PGADrawPolyline([Q1,Q3],'r','LineWidth',2);");
     PGADrawPolyline([P1,P3],'r','LineWidth',2); PGADrawPolyline([Q1,Q3],'r','LineWidth',2);
     disp(">> PGADrawPolyline([P1,Q2],'k','LineWidth',2); PGADrawPolyline([P1,Q3],'k','LineWidth',2);");
     PGADrawPolyline([P1,Q2],'k','LineWidth',2); PGADrawPolyline([P1,Q3],'k','LineWidth',2);
     disp(">> PGADrawPolyline([P2,Q1],'k','LineWidth',2); PGADrawPolyline([P2,Q3],'k','LineWidth',2);");
     PGADrawPolyline([P2,Q1],'k','LineWidth',2); PGADrawPolyline([P2,Q3],'k','LineWidth',2);
     disp(">> PGADrawPolyline([P3,Q1],'k','LineWidth',2); PGADrawPolyline([P3,Q2],'k','LineWidth',2);");
     PGADrawPolyline([P3,Q1],'k','LineWidth',2); PGADrawPolyline([P3,Q2],'k','LineWidth',2);

     GAprompt;

     disp('>> H3 = IntersectLines(join(P1,Q2),join(P2,Q1));');
     H3 = IntersectLines(join(P1,Q2),join(P2,Q1));
     disp('>> H2 = IntersectLines(join(P1,Q3),join(P3,Q1));');
     H2 = IntersectLines(join(P1,Q3),join(P3,Q1));
     disp('>> H1 = IntersectLines(join(P2,Q3),join(P3,Q2));');
     H1 = IntersectLines(join(P2,Q3),join(P3,Q2));
     disp(">> draw(H1,'g'); draw(H2,'g'); draw(H3,'g')");
     draw(H1,'g'); draw(H2,'g'); draw(H3,'g')
     disp(">> PGADrawPolyline([H1,H3],'b','LineWidth',3);");
     PGADrawPolyline([H1,H3],'b','LineWidth',3);
     disp('>> view([0 90])');
     view([0 90])
     GAprompt;
  disp(' ');    disp('End of PGAblock sequence.  Returning to Matlab.');
elseif ( GAn == 3 )
     %PGAaddendum(3)
     disp(">> %MORLEY'S TRIANGLE");
     %MORLEY'S TRIANGLE
     disp(">> clf");
     clf
     disp(">> P1 = gapoint(0,0,0); P2 = gapoint(1,0,0); P3 = gapoint(0,1,0);");
     P1 = gapoint(0,0,0); P2 = gapoint(1,0,0); P3 = gapoint(0,1,0);
     disp(">> draw(P1,'r'); draw(P2,'g'); draw(P3,'b');");
     draw(P1,'r'); draw(P2,'g'); draw(P3,'b');
     disp(">> PGADrawPolyline([P1,P2,P3,P1], 'k','LineWidth',2);");
     PGADrawPolyline([P1,P2,P3,P1], 'k','LineWidth',2);
     disp(">> GAview([0,90]);");
     GAview([0,90]);
     disp(">> J12=join(P1,P2); J12 = J12/norm(J12);");
     J12=join(P1,P2); J12 = J12/norm(J12);
     disp(">> J23=join(P2,P3); J23 = J23/norm(J23);");
     J23=join(P2,P3); J23 = J23/norm(J23);
     disp(">> J31=join(P3,P1); J31 = J31/norm(J31);");
     J31=join(P3,P1); J31 = J31/norm(J31);
     disp(">> R = -J12/J31;");
     R = -J12/J31;
     disp(">> R1 = gexp(GAZ(glog(R))/6); R1i = inverse(R1);");
     R1 = gexp(GAZ(glog(R))/6); R1i = inverse(R1);
     disp(">> R = -J23/J12; Ri=inverse(R);");
     R = -J23/J12; Ri=inverse(R);
     disp(">> R2 = gexp(GAZ(glog(R))/6); R2i = inverse(R2);");
     R2 = gexp(GAZ(glog(R))/6); R2i = inverse(R2);
     disp(">> R = -J31/J23;");
     R = -J31/J23;
     disp(">> R3 = gexp(GAZ(glog(R))/6); R3i = inverse(R3);");
     R3 = gexp(GAZ(glog(R))/6); R3i = inverse(R3);

     GAprompt;

     disp(">> L31 = R3i*J31*R3; draw(L31,'c')");
     L31 = R3i*J31*R3; draw(L31,'c')
     disp(">> L13 = R1*J31*R1i; draw(L13,'c');");
     L13 = R1*J31*R1i; draw(L13,'c');
     disp(">> EP1 = IntersectLines(L13,L31);");
     EP1 = IntersectLines(L13,L31);

     GAprompt;

     disp(">> L23 = R2i*J23*R2; draw(L23,'m')");
     L23 = R2i*J23*R2; draw(L23,'m')
     disp(">> L32 = R3*J23*R3i; draw(L32,'m')");
     L32 = R3*J23*R3i; draw(L32,'m')
     disp(">> EP2 = IntersectLines(L23,L32);");
     EP2 = IntersectLines(L23,L32);
     disp(">> L12 = R1i*J12*R1; draw(L12,'y');");
     L12 = R1i*J12*R1; draw(L12,'y');
     disp(">> L21 = R2*J12*R2i; draw(L21,'y');");
     L21 = R2*J12*R2i; draw(L21,'y');
     disp(">> EP3 = IntersectLines(L12,L21);");
     EP3 = IntersectLines(L12,L21);
     disp(">> PGADrawPolyline([EP1,EP2,EP3,EP1],'k','LineWidth',3)");
     PGADrawPolyline([EP1,EP2,EP3,EP1],'k','LineWidth',3)
     disp(">> view([0 90])");
     view([0 90])

     GAprompt;
     disp(' ');    disp('End of PGAblock sequence.  Returning to Matlab.');
elseif ( GAn == 4 ) 
  GAps = 'PGAblock >> ';
     %PGAblock(4)
     disp('        % REFLECTION IN A PLANE');
     % REFLECTION IN A PLANE
     disp('        n1 = e1-0.5*e0;    % construct a plane to reflect through');
     n1 = e1-0.5*e0;
     disp('        n2 = (e1+e2+e3)/norm(e1+e2+e3)-0.5*e0; % construct a plane to reflect');
     n2 = (e1+e2+e3)/norm(e1+e2+e3)-0.5*e0;
     disp('        Pt1 = n2^e1^e2;                 % construct point on plane n2');
     Pt1 = n2^e1^e2;                 % construct point on plane n2
     disp('        Pt1 = -Pt1/inner(Pt1,e1^e2^e3); % normalize the point');
     Pt1 = -Pt1/inner(Pt1,e1^e2^e3); % normalize the point
     disp('        Pt1r = -n1*Pt1*n1;               % compute its reflection');
      
     Pt1r = -n1*Pt1*n1;
     disp("        pclf; draw(n1); draw(Pt1); draw(Pt1r,'r'); GAview([-15 21]);%%");
     pclf; draw(n1); draw(Pt1); draw(Pt1r,'r'); GAview([-15 21]);%%
     GAprompt;
     disp('        n2r =  -n1*n2*n1;               % reflect plane n2 through plane n1');
     n2r = -n1*n2*n1;
     disp("        draw(n2,'b'); draw(n2r,'m'); GAview([-5,50]);");
     draw(n2,'b'); draw(n2r,'m'); GAview([-5,50]);    %%
     GAprompt;
     disp('        L = e1^e2;                      % construct a line');
     L = e1^e2;  % Construct a line
     disp('        Lr = n1*L*n1;                  % reflect line L through plane n1');
     Lr = n1*L*n1;
     disp("        draw(L); draw(Lr); GAview([-5,50]);   %%");
     draw(L); draw(Lr); GAview([-5,50]);
     GAprompt;
     disp(' ');    disp('End of GAblock sequence.  Returning to Matlab.');
elseif ( GAn == 5 ) 
  GAps = 'PGAblock >> ';
     disp("         % TRANSLATION");
     % TRANSLATION");
     disp("         n1=(e1+e2+e3)/norm(e1+e2+e3)-0.5*e0;");
     n1=(e1+e2+e3)/norm(e1+e2+e3)-0.5*e0;
     disp("         n2=(e1+e2+e3)/norm(e1+e2+e3)-1.5*e0;");
     n2=(e1+e2+e3)/norm(e1+e2+e3)-1.5*e0;
     disp("         Pt = (1-e0*(e1+e3))*e123;");
     Pt = (1-e0*(e1+e3))*e123;
     disp("         pclf; draw(n1); draw(n2,'b'); draw(Pt);  GAview([30,36]) %%");
     pclf; draw(n1); draw(n2,'b'); draw(Pt);  GAview([30,36]) %% 
     GAprompt;
     disp("         Ptr1 = n1*Pt*n1; draw(Ptr1,'m');  GAview([30,36])       %% First reflection ");
     Ptr1 = n1*Pt*n1; draw(Ptr1,'m');  GAview([30,36])       %% First reflection
     GAprompt
     disp("         Ptr2 = n2*Ptr1*n2; draw(Ptr2,'r'); GAview([30,36])      %% Second reflection ");
     Ptr2 = n2*Ptr1*n2; draw(Ptr2,'r'); GAview([30,36])      %% Second reflection
     disp(' ');
     GAps = 'Refer to the tutorial before continuing >> ';
     GAprompt
     disp("         Tt = n2*n1");
     Tt = n2*n1
     disp("         t  = 2*(1.5-0.5)*(e1+e2+e3)/norm(e1+e2+e3);");
     t  = 2*(1.5-0.5)*(e1+e2+e3)/norm(e1+e2+e3);
     disp("         Tt = 1-e0*t/2                                           %%  Observe the two Tt's are the same");
     Tt = 1-e0*t/2                                           %%  Observe the two Tt's are the same
     disp(' ');
     GAps = 'Refer to the tutorial before continuing >> ';
     GAprompt;
     disp("         L = e1^e2; draw(L); GAview([30,36]);                    %% A line");
     L = e1^e2; draw(L); GAview([30,36]);                    %% A line
     GAprompt
     disp("         Lt = Tt*L*inverse(Tt); draw(Lt); GAview([30,36]);       %% The line translated    ");
     Lt = Tt*L*inverse(Tt); draw(Lt); GAview([30,36]);       %% The line translated    
     GAprompt
  disp(' ');    disp('End of GAblock sequence.  Returning to Matlab.');
elseif ( GAn == 6 ) 
  GAps = 'GAblock >> ';
     disp("         % ROTATIONS");
     % ROTATIONS
     disp("         n1=e1-0.25*e0;");
     n1=e1-0.25*e0;
     disp("         n2=(5*e1+e2+e3)/norm(5*e1+e2+e3)-0.25*e0;");
     n2=(5*e1+e2+e3)/norm(5*e1+e2+e3)-0.25*e0;
     disp("         Tr = n2*n1;");
     Tr = n2*n1;
     disp("         Pt = (1-e0*(e1+e3))*e123;");
     Pt = (1-e0*(e1+e3))*e123;
     disp("         pclf; draw(n1); draw(n2,'b'); draw(Pt);");
     pclf; draw(n1); draw(n2,'b'); draw(Pt);
     GAprompt
     disp("         for i=1:10");
     for i=1:10
       fprintf(1,"           i=%d\n",i);
       disp("           Pt = zeroepsilons(Tr*Pt*inverse(Tr));");
       Pt = zeroepsilons(Tr*Pt*inverse(Tr));
       disp("           draw(Pt,'r');");
       draw(Pt,'r');
       pause(1)
     end
     disp("         end");
     disp("         %%");
     %%
     GAprompt;

  disp(' ');    disp('End of GAblock sequence.  Returning to Matlab.');
elseif ( GAn == 7 ) 
  GAps = 'GAblock >> ';
     disp("     % SCREWS");
     clf
     % SCREWS
     disp("     tv = (e1+e2+0.5*e3)/10;");
     tv = (e1+e2+0.5*e3)/10;
     disp("     P =gapoint(1,0,1);");
     P =gapoint(1,0,1);
     disp("     L = P.*tv/norm(tv);");
     L = P.*tv/norm(tv);
     disp("     draw(P); draw(e123); draw(L);");
     draw(P); draw(e123); draw(L);

     disp(" ");
     disp("     Q =gapoint(0,1,0);");
     Q =gapoint(0,1,0);
     disp("     draw(Q,'b');");
     draw(Q,'b');

     disp(" ");
     disp("     T = gexp(-e0*tv/2);");
     T = gexp(-e0*tv/2);
     disp("     Ti = gexp(e0*tv/2);");
     Ti = gexp(e0*tv/2);
     GAprompt;

     disp("for i=1:5");
     for i=1:5
         fprintf(1,"           i=%d\n",i);
         disp("    	Q = T*Q*T");
         Q = T*Q*Ti;
         disp("    	draw(Q,'r');");
         draw(Q,'r');
	 pause(0.2);
     end;
     disp("end");
     GAprompt;

     disp(" ");
     disp("    Q =gapoint(0,1,0);");
     Q =gapoint(0,1,0);
     disp("    RL = gexp(-pi/12*L/2);");
     RL = gexp(-pi/12*L/2);
     disp("    RLi = gexp(pi/12*L/2);");
     RLi = gexp(pi/12*L/2);
     GAprompt;

     disp("for i=1:1:23");
     for i=1:23
         if i<3
            fprintf(1,"           i=%d\n",i);
            disp("	Q = RL*Q*RLi;");
	 end
	 Q = RL*Q*RLi;
	 if i<3
            disp("	draw(Q,'g');");
	 end
	 if i>=3 
	    fprintf(1,".");
	 end
	 draw(Q,'g');
	 pause(0.1);
     end
     disp(" ");
     disp("end");
     GAprompt;

     disp(" ");
     disp("R =gapoint(0,0,0);");
     R =gapoint(0,0,0);
     disp("Sc = T*RL; Sci = RLi*Ti;");
     Sc = T*RL; Sci = RLi*Ti;
     disp("draw(R,'y');");
     draw(R,'y');
     GAprompt;

     disp("for i=1:30");
     for i=1:30
         if i<3
           fprintf(1,"           i=%d\n",i);
           disp("	R = GAZ(Sc*R*Sc);");
         end
	 R = GAZ(Sc*R*Sci);
	 if i<3
            disp("	draw(R,'m');");
	 end
	 if i>=3
	    fprintf(1,".");
	 end
	 draw(R,'m');
	 pause(0.05);
     end
     fprintf(1,"\n");
     disp("end");
     GAprompt;
     disp("    GAview([-45,-19])");
     GAview([-45,-19])
     GAprompt;
     disp("    GAorbiter(380,6)");
     GAorbiter(380,6)
end
catch ; end

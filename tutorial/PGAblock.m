function PGAblock(GAn)
%GAblock: run sample code in tutorial.

% PGABLE, Copyright (c) 2024, University of Waterloo
% Copying, use and development for non-commercial purposes permitted.
%          All rights for commercial use reserved; for more information
%          contact Stephen Mann (smann@uwaterloo.ca)
%
%          This software is unsupported.

try
if ( GAn == 1 ) 
  GAps = 'PGAblock >> ';
     disp('>>      % POINTS');
     disp('>>      n1 = e3-.5*e0;');
     n1 = e3-.5*e0;
     disp('>>      n2 = e1-0.25*e0;');
     n2 = e1-0.25*e0;
     disp('>>      n3 = e2+0.75*e0;');
     n3 = e2+0.75*e0;
     disp('>>      clf; draw(n1); draw(n2,’r’); draw(n3,’b’);  %%');
     clf; draw(n1); draw(n2,'r'); draw(n3,'b');
     GAprompt;
     disp('>>      P = n1^n2^n3');
     P = n1^n2^n3
     disp('>>      draw(P)');
     draw(P)
     disp(' ');    disp('End of PGAblock sequence.  Returning to Matlab.');
elseif ( GAn == 2 ) 
  GAps = 'PGAblock >> ';
     disp('        % INTERSECTING NON-ORTHOGONAL PLANES');
     % INTERSECTING NON-ORTHOGONAL PLANES
     disp('        n1 = (e1+e2)/norm(e1+e2) -0.25*e0;');
     n1 = (e1+e2)/norm(e1+e2) -0.25*e0;
     disp('        n2 = e2 +0.5*e0;');
     n2 = e2 +0.5*e0;
     disp('        n3 = e3 -0.25*e0;');
     n3 = e3 -0.25*e0;
     disp("        clf; draw(n1); draw(n2,'r'); draw(n3,'b');");
     clf; draw(n1); draw(n2,'r'); draw(n3,'b');
     disp('        P = n1^n2^n3');
     P = n1^n2^n3
     disp('        draw(P)');
     draw(P)
     disp('        GAview([65 30])  %%');
     GAview([65 30])
     GAprompt;
  disp(' ');    disp('End of PGAblock sequence.  Returning to Matlab.');
elseif ( GAn == 3 )
     disp("         % POINTS ON PLANES");
     % POINTS ON PLANES
     disp("         n = e3-0.5*e0;");
     n = e3-0.5*e0;
     disp("         clf; draw(n); ");
     clf; draw(n); 
     disp("         x1=0.5*e3;    P1 = (1-e0*x1)*e123;");
     x1=0.5*e3;    P1 = (1-e0*x1)*e123;
     disp("         x2=e1+0.5*e3; P2 = (1-e0*x2)*e123;");
     x2=e1+0.5*e3; P2 = (1-e0*x2)*e123;
     disp("         draw(P1); draw(P2); view([-58, 12]); %%");
     draw(P1); draw(P2); view([-58, 12]); %%
     GAprompt;
     disp("         x3=e3;        P3 = (1-e0*x3)*e123;");
     x3=e3;        P3 = (1-e0*x3)*e123;
     disp("         x4=0;         P4 = (1-e0*x4)*e123;");
     x4=0;         P4 = (1-e0*x4)*e123;
     disp("         draw(P3); draw(P4); view([-58, 12]); %%");
     draw(P3); draw(P4); view([-58, 12]); %%
     GAprompt;
  disp(' ');    disp('End of PGAblock sequence.  Returning to Matlab.');
elseif ( GAn == 4 ) 
  GAps = 'PGAblock >> ';
     disp('        %PGAblock(3)');
     %PGAblock(3)
     disp('        % REFLECTION IN A PLANE');
     % REFLECTION IN A PLANE
     disp('        n1 = e1-0.5*e0;');
     n1 = e1-0.5*e0;
     disp('        n2 = (e1+e2+e3)/norm(e1+e2+e3)-0.5*e0;');
     n2 = (e1+e2+e3)/norm(e1+e2+e3)-0.5*e0;
     disp('        Pt1 = n2^e1^e2;                 % construct point on plane n2');
     Pt1 = n2^e1^e2;                 % construct point on plane n2
     disp('        Pt1 = -Pt1/inner(Pt1,e1^e2^e3); % normalize the point');
     Pt1 = -Pt1/inner(Pt1,e1^e2^e3); % normalize the point
     disp('        Pt1r = n1*Pt1*n1;');
     Pt1r = n1*Pt1*n1;
     disp("        clf; draw(n1); draw(Pt1); draw(Pt1r,'r'); %%");
     clf; draw(n1); draw(Pt1); draw(Pt1r,'r'); %%
     GAprompt;
     disp('        n2r =  -n1*n2*n1;');
     n2r = -n1*n2*n1;
     disp("        draw(n2,'b'); draw(n2r,'m'); GAview([-5,50]);");
     draw(n2,'b'); draw(n2r,'m'); GAview([-5,50]);    %%
     GAprompt;
     disp('        L = e1^e2;  % Construct a line');
     L = e1^e2;  % Construct a line
     disp('        Lr = -n1*L*n1;');
     Lr = -n1*L*n1;
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
     disp("         clf; draw(n1); draw(n2,'b'); draw(Pt);  GAview([30,36]) %%");
     clf; draw(n1); draw(n2,'b'); draw(Pt);  GAview([30,36]) %% 
     GAprompt;
     disp("         Ptr1 = n1*Pt*n1; draw(Ptr1,'m');  GAview([30,36])       %% First reflection ");
     Ptr1 = n1*Pt*n1; draw(Ptr1,'m');  GAview([30,36])       %% First reflection
     GAprompt
     disp("         Ptr2 = n2*Ptr1*n2; draw(Ptr2,'r'); GAview([30,36])      %% Second reflection ");
     Ptr2 = n2*Ptr1*n2; draw(Ptr2,'r'); GAview([30,36])      %% Second reflection
     GAprompt
     disp("         Tt = n2*n1");
     Tt = n2*n1
     disp("         t  = 2*(1.5-0.5)*(e1+e2+e3)/norm(e1+e2+e3);");
     t  = 2*(1.5-0.5)*(e1+e2+e3)/norm(e1+e2+e3);
     disp("         Tt = 1-e0*t/2                                           %%  Observe the two Tt's are the same");
     Tt = 1-e0*t/2                                           %%  Observe the two Tt's are the same
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
     disp("         clf; draw(n1); draw(n2,'b'); draw(Pt);");
     clf; draw(n1); draw(n2,'b'); draw(Pt);
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
end
catch ; end

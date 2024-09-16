classdef OGA < GA
    % OGA  is a child class of GA for elements of Projective/Plane-based Geometric Algebra.
    %
    % See also GA, PGA, CGA.
    properties (Access = private)
        % A 1x8 matrix of real numbers corresponding to the coefficients of entries 1, e1, e2, e3, e12, e13, e23, e123. 
        m
    end

    %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%
    %           Settings           %
    %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%

    methods (Static = true)
        % TODO: Consider removing this functionality, perhaps?
        % TODO: Consider renaming to just 'signature'
        function [S1, S2, S3] = signature(sign1, sign2, sign3)
            persistent signature1;
            persistent signature2;
            persistent signature3;

            if isempty(signature1)
                signature1 = 1;
                signature2 = 1;
                signature3 = 1;
            end

            if nargin == 3
               signature1 = sign1;
               signature2 = sign2;
               signature3 = sign3;
            end

            S1 = signature1;
            S2 = signature2;
            S3 = signature3;
        end
    end

    %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%
    %      Protected methods       %
    %         (non-static)         %
    %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%

    methods (Access = protected)
        function b = GAisa_(A, t)
            nzm = A.m ~= 0;
            if strcmp(t, 'double') || strcmp(t, 'scalar') 
                b = sum(nzm(2:8)) == 0;
            elseif strcmp(t,'vector') 
                b = sum(nzm(5:8)) == 0 & nzm(1) == 0;
            elseif strcmp(t,'bivector') 
                b = sum(nzm(1:4)) == 0 & nzm(8) == 0;
            elseif strcmp(t,'trivector') 
                b = sum(nzm(1:7)) == 0;
            elseif strcmp(t,'multivector')
                b = sum( [sum(nzm(1)) sum(nzm(2:4)) sum(nzm(5:7)) nzm(8)] ~= 0) > 1;
            else
                b = false;
            end 
        end

        function r = double_(A)
            if GAisa_(A, 'scalar')
                r = A.m(1);
            else
                error('Can only convert a scalar OGA object to a double.');
            end
        end

        function R = plus_(A, B)
            R = OGA(A.m + B.m);
        end

        function R = minus_(A, B)
            R = OGA(A.m - B.m);
        end

        function R = uminus_(A)
            R = OGA(-A.m);
        end

        function mr = productleftexpand_(A)
            [S1, S2, S3] = OGA.signature();
            S12 = S1*S2;
            S13 = S1*S3;
            S23 = S2*S3;
            S123 = S12*S3;

            scal = A.m(1);
            E1 = A.m(2);
            E2 = A.m(3);
            E3 = A.m(4);
            E12 = A.m(5);
            E13 = A.m(6);
            E23 = A.m(7);
            E123 = A.m(8);

            mr = [scal     S1*E1      S2*E2      S3*E3     -S12*E12   -S13*E13   -S23*E23  -S123*E123 ;
                 E1       scal       S2*E12     S3*E13    -S2*E2     -S3*E3     -S23*E123  -S23*E23  ;
                 E2      -S1*E12     scal       S3*E23     S1*E1      S13*E123  -S3*E3      S13*E13  ;
                 E3      -S1*E13    -S2*E23     scal      -S12*E123   S1*E1      S2*E2     -S12*E12  ;
                 E12     -E2         E1         S3*E123    scal       S3*E23    -S3*E13    S3*E3     ;
                 E13     -E3        -S2*E123    E1        -S2*E23     scal       S2*E12   -S2*E2     ;
                 E23      S1*E123   -E3         E2         S1*E13    -S1*E12     scal      S1*E1     ;
                 E123     E23       -E13        E12        E3        -E2         E1        scal      ];
        end

        function R = product_(A, B)
            R = OGA(productleftexpand_(A)*B.m);
        end 

        function mr = outerleftexpand_(A)
            scal = A.m(1);
            E1 = A.m(2);
            E2 = A.m(3);
            E3 = A.m(4);
            E12 = A.m(5);
            E13 = A.m(6);
            E23 = A.m(7);
            E123 = A.m(8);

            mr = [scal     0       0       0       0      0       0       0   ;
                 E1      scal     0       0       0      0       0       0   ;
                 E2       0      scal     0       0      0       0       0   ;
                 E3       0       0      scal     0      0       0       0   ;
                 E12     -E2      E1      0      scal    0       0       0   ;
                 E13     -E3      0       E1      0     scal     0       0   ;
                 E23      0      -E3      E2      0      0      scal     0   ;
                 E123    E23     -E13    E12      E3    -E2      E1     scal];
        end

        function R = outer_(A, B)
            R = OGA(outerleftexpand_(A)*B.m);
        end

        function R = inner_(A, B)
            [S1, S2, S3] = OGA.signature();
            S12 = S1*S2;
            S13 = S1*S3;
            S23 = S2*S3;
            S123 = S12*S3;

            [scal, E1, E2, E3, E12, E13, E23, E123] = decompose_(A);

            C1 = S1*E1;
            C2 = S2*E2;
            C3 = S3*E3;
            C12 = S12*E12;
            C13 = S13*E13;
            C23 = S23*E23;
            C123 = S123*E123;

            % TODO: write code to generate this matrix.
            M = [scal   C1    C2    C3  -C12  -C13  -C23  -C123 ;

                 C1    scal   0     0    -C2   -C3    0   -C23  ;
                 C2   -C12   scal   0     C1    0    -C3   C13  ;
                 C3   -C13  -C23   scal   0     C1    C2  -C12  ;

                 C12    0     0    C123  scal   0     0     C3  ;
                 C13    0   -C123   0     0    scal   0    -C2  ;
                 C23   C123   0     0     0     0    scal   C1  ;

                 C123   0     0     0     0     0     0    scal ];
            R = OGA(M*B.m);
        end

        function [scal, E1, E2, E3, E12, E13, E23, E123] = decompose_(A)
            scal = A.m(1);
            E1 = A.m(2);
            E2 = A.m(3);
            E3 = A.m(4);
            E12 = A.m(5);
            E13 = A.m(6);
            E23 = A.m(7);
            E123 = A.m(8);
        end

        function R = leftcontraction_(A, B)
            % TODO: Implement.
            
            error('Not yet implemented.')
        end

        function R = rightcontraction_(A, B)
            % TODO: Implement.

            error('Not yet implemented.')
        end

        function R = inverse_(A)
            rm = productleftexpand_(A);
            if rcond(rm) <= eps
                error('Inverse of %s does not exist.', char(A))
            end

           R = OGA(rm\OGA(1).m);
        end

        function R = divide_(A, B)
            R = A * inverse_(B);
        end

        function r = norm_(A)
            B = double_(grade_(A.product_(reverse_(A)), 0));
            if B > 0
                r = sqrt(B);
            else
                r = sqrt(-B);
            end
        end

        function r = vnorm_(A)
            error('vnorm does not exist in OGA. Try norm.')
        end

        function R = normalize_(A)
            R = A / norm_(A);
        end

        function b = eq_(A, B)
            b = norm_(A - B) < GA.epsilon_tolerance;
        end

        function b = eeq_(A, B)
            b = all(A.m == B.m);
        end

        function b = ne_(A, B)
            b = ~eq_(A, B);
        end

        function R = dual_(A)
            R = A/I3(OGA);
        end

        function R = inversedual_(A)
            R = A*I3(OGA);
        end

        function R = hodgedual_(A)
            error('Hodge dual cannot be performed on OGA elements.');
        end

        function R = inversehodgedual_(A)
            error('Inverse Hodge dual cannot be performed on OGA elements.');
        end

        function R = jmap_(A)
            error('Jmap cannot be performed on OGA elements.');
        end

        function R = reverse_(A)
            R = OGA([A.m(1);A.m(2);A.m(3);A.m(4);-A.m(5);-A.m(6);-A.m(7);-A.m(8)]);
        end

        function R = zeroepsilons_(A)
            R = A;
            for i = 1:8
                if abs(R.m(i)) < GA.epsilon_tolerance
                    R.m(i) = 0;
                end
            end
        end

        function R = join_(A, B)
            % TODO: Implement.

            error('Not yet implemented.');
        end

        function R = meet_(A, B)
            % TODO: Implement.

            error('Not yet implemented.');
        end

        function R = conjugate_(A)
            R = OGA([A.m(1); -A.m(2); -A.m(3); -A.m(4); -A.m(5); -A.m(6);  -A.m(7);  A.m(8)]);
        end

        function R = gradeinvolution_(A)
            R = OGA([A.m(1); -A.m(2); -A.m(3); -A.m(4); A.m(5); A.m(6); A.m(7); -A.m(8)]);
        end

        function R = grade_(A, n)
            % grade(A, n) return the part of an object of a particular grade.
            %  Return the part of A of grade n.
            %  If n is omitted, return the grade of A (-1 if A is of mixed grade) 
            
            if nargin == 1
                if A.m(1) ~= 0
                    if sum(abs(A.m(2:8))) == 0
                        R = 0;
                    else
                        R = -1;
                    end
                elseif sum(abs(A.m(2:4))) ~= 0
                    if sum(abs(A.m(5:8))) == 0
                        R = 1;
                    else
                        R = -1;
                    end
                elseif sum(abs(A.m(5:7))) ~= 0
                    if A.m(8) == 0
                        R = 2;
                    else
                        R = -1;
                    end
                elseif A.m(8) ~= 0
                    R = 3;
                else
                    R = -1;
                end
            else
                if n == 0
                    R = OGA([A.m(1);0;0;0;0;0;0;0]);
                elseif n == 1
                    R = OGA([0;A.m(2);A.m(3);A.m(4);0;0;0;0]);
                elseif n == 2
                    R = OGA([0;0;0;0;A.m(5);A.m(6);A.m(7);0]);
                elseif n == 3
                    R = OGA([0;0;0;0;0;0;0;A.m(8)]);
                else
                    R = OGA(0);
                end
            end
        end

        function b = isgrade_(A, g)
            b = false;

            if g == 0
                if sum(abs(A.m(2:8))) == 0
                    b = true;
                end
            elseif g == 1
                if sum(abs([A.m(1);A.m(5:8)])) == 0
                    b = true;
                end
            elseif g == 2
                if sum(abs([A.m(1:4);A.m(8)])) == 0
                    b = true;
                end
            elseif g == 3
                if sum(abs(A.m(1:7))) == 0
                    b = true;
                end
            elseif g == -1
                z = A.m == 0;
                z0 = z(1);
                z1 = sum(z(2:4)) ~= 0;
                z2 = sum(z(5:7)) ~= 0;
                z3 = z(8);
                % TODO: What does the comment below mean?
                % Note that the test treats 0 is a multivector!
                if z0 + z1 + z2 + z3 ~= 1
                    b = true;
                end
            else
                error('isgrade: invalid grade.');
            end
        end

        % function r = unit(A)
        %     if isGrade(A,0)
        %         r = unit(A.m(1));
        %     elseif isGrade(A,1)
        %         s = sqrt(A.m(2)*A.m(2)+A.m(3)*A.m(3)+A.m(4)*A.m(4));
        %         r = A/s;
        %     elseif isGrade(A,2)
        %         s = sqrt(A.m(5)*A.m(5)+A.m(6)*A.m(6)+A.m(7)*A.m(7));
        %         r = A/s;
        %     elseif isGrade(A,3)
        %         r = A/abs(A.m(8));
        %     else
        %         error('Unit can only be applied to blades.');
        %     end
        % end

        % function r = blade(A)
        %     % blade(A) : return a blade made from the largest portion of a multivector.
        %     A = OGA(A);

        %     s(1) = abs(A.m(1));
        %     s(2) = sqrt(sum(abs(A.m(2:4))));
        %     s(3) = sqrt(sum(abs(A.m(5:7))));
        %     s(4) = abs(A.m(8));
        %     if s(1)>s(2) && s(1)>s(3) && s(1)>s(4)
        %       r = OGA.returnOGA(A.m(1));
        %     elseif s(2)>s(3) && s(2)>s(4)
        %       r = OGA.returnOGA([0; A.m(2); A.m(3); A.m(4); 0; 0; 0; 0]);
        %     elseif s(3)>s(4)
        %       r = OGA.returnOGA([0; 0; 0; 0; A.m(5); A.m(6); A.m(7); 0]);
        %     else
        %       r = OGA.returnOGA([0; 0; 0; 0; 0; 0; 0; A.m(8)]);
        %     end
        % end

        function R = gexp_(A)
            %gexp(A): Computes the geometric product exponential of a multivector.
            rm = productleftexpand_(A);
            E = expm(rm);
            R = OGA(E(1:8,1));
        end

        % function r = wexp_(A)
        %     %GAwexp(A): Gives the wedge product exponential of a GA objects.
        
        %     M =[A.m(1)    0       0       0       0       0       0       0;
        %     A.m(2)  A.m(1)    0       0       0       0       0       0;
        %     A.m(3)    0     A.m(1)    0       0       0       0       0;
        %     A.m(4)    0       0     A.m(1)    0       0       0       0;
        %     A.m(5) -A.m(3)  A.m(2)    0     A.m(1)    0       0       0;
        %     A.m(6)    0    -A.m(4)  A.m(3)    0     A.m(1)    0       0;
        %     A.m(7)  A.m(4)    0    -A.m(2)    0       0     A.m(1)    0;
        %     A.m(8)  A.m(6)  A.m(7)  A.m(5)  A.m(4)  A.m(2)  A.m(3)  A.m(1)];
        %     E = expm(M);
        %     r = OGA(E(1:8,1));
        % end

        function R = glog_(A)
            rm = productleftexpand_(A);
            L = logm(rm);
            R = OGA(L(1:8, 1));
        end

        function R = sqrt_(A)
            rm = productleftexpand_(A);
            S = sqrtm(rm);
            R = PGA(S(1:8, 1));
            % TODO: fix this. Gets i on -1.
        end

        function r = getx_(A)
            % Position of e1
            r = A.m(2);
        end
        
        function r = gety_(A)
            % Position of e2
            r = A.m(3);
        end
        
        function r = getz_(A)
            % Position of e3
            r = A.m(4);
        end

        function s = char_(p)
            if ~any(p.m(:))
                s = '0';
                return;
            end

            pl = '';
            s = '';
            if p.m(1) ~= 0
                s = [s pl num2str(p.m(1))];
                pl = ' + ';
            end

            [s, pl] = GA.charify_val_(p.m(2), 'e1', s, pl);
            [s, pl] = GA.charify_val_(p.m(3), 'e2', s, pl);
            [s, pl] = GA.charify_val_(p.m(4), 'e3', s, pl);
                        
            if GA.compact_notation()
                if ~GA.increasing_order()
                    [s, pl] = GA.charify_val_(p.m(7), 'e23', s, pl);
                    [s, pl] = GA.charify_val_(-p.m(6), 'e31', s, pl);
                    [s, pl] = GA.charify_val_(p.m(5), 'e12', s, pl);
                else
                    [s, pl] = GA.charify_val_(p.m(5), 'e12', s, pl);
                    [s, pl] = GA.charify_val_(p.m(6), 'e13', s, pl);
                    [s, pl] = GA.charify_val_(p.m(7), 'e23', s, pl);
                end

                if GA.compact_pseudoscalar()
                    [s, pl] = GA.charify_val_(p.m(8), 'I3', s, pl);
                else
                    [s, pl] = GA.charify_val_(p.m(8), 'e123', s, pl);
                end
            else
                if ~GA.increasing_order()
                    [s, pl] = GA.charify_val_(p.m(5), 'e1^e2', s, pl);
                    [s, pl] = GA.charify_val_(p.m(6), 'e1^e3', s, pl);
                    [s, pl] = GA.charify_val_(p.m(7), 'e2^e3', s, pl);
                else 
                    [s, pl] = GA.charify_val_(p.m(7), 'e2^e3', s, pl);
                    [s, pl] = GA.charify_val_(-p.m(6), 'e3^e1', s, pl);
                    [s, pl] = GA.charify_val_(p.m(5), 'e1^e2', s, pl);
                end

                if GA.compact_pseudoscalar()
                    [s, pl] = GA.charify_val_(p.m(8), 'I3', s, pl);
                else
                    [s, pl] = GA.charify_val_(p.m(8), 'e1^e2^e3', s, pl);
                end
            end
            
            if strcmp(pl, ' ')
                s = '     0';
            end
        end
    end
    
    methods (Access = public)
        function obj = OGA(m0, m1, m2, m3)
            % Constructor.
        
            % TODO: Leaved detailed description.

            if nargin == 0
                obj = OGA(0);
            elseif nargin == 4
                if m1 == 0
                    m1 = zeros(3, 1);
                end
                if m2 == 0
                    m2 = zeros(3, 1);
                end
                obj.m = [m0; 
                         m1(1); m1(2); m1(3)
                         m2(1); m2(2); m2(3);  
                         m3];
            elseif nargin == 1
                if isa(m0, 'OGA')
                    obj = m0;
                elseif isa(m0, 'double')
                    if size(m0, 1) == 1 & size(m0, 2) == 1
                        % User has provided a scalar
                        obj.m = [m0; 
                                 0; 0; 0; 
                                 0; 0; 0; 
                                 0];
                    elseif size(m0, 1) == 1 & size(m0, 2) == 8
                        % User has provided a column vector
                        obj.m = m0';
                    elseif size(m0, 1) == 8 & size(m0, 2) == 1
                        % User has prodived a row vector
                        obj.m = m0;
                    else
                        error('Bad OGA argument: Unexpected array size.\nExpected size is either 1x1, 8x1 or 1x8.\nCurrent size is: %dx%d', size(m0, 1), size(m0, 2))
                    end
                else
                    error('Bad OGA argument: Invalid input type. Class is currently: %s.', class(m0))
                end
            else 
                error('Bad OGA argument: Invalid number of arguments.')
            end
        end

        % Returns the matrix for the OGA object. For debugging purposes.
        % TODO: add to change notes: method m in GA is now called matrix in OGA.
        function rm = matrix(self)
            rm = self.m;
        end

        function b = GAisa(A, t)
            arguments
                A OGA;
                t (1, 1) string;
            end
            
            b = GAisa_(A, t);
        end

        function R = OGAcast(A)
            R = A;
        end

        function R = PGAcast(A)
            scal = A.m(1);
            E1   = A.m(2);
            E2   = A.m(3);
            E3   = A.m(4);
            E12  = A.m(5);
            E13  = A.m(6);
            E23  = A.m(7);
            E123 = A.m(8);
            R = PGA(scal, [0, E1, E2, E3], [0, 0, 0, E12, E13, E23], [0, 0, 0, E123], 0);
        end

        function draw(A, c)
            arguments
                A OGA;
                c = [];
            end

            A = zeroepsilons_(A);
            
            if GAisa(A, 'vector')
                if isempty(c)
                    c = 'b';
                end
                % TODO: make a proper way of converting a vector to a point
                h = GAScene.drawarrow(origin(PGA), ihd(PGAcast(A) + e0(PGA)), c);
                GAScene.additem(GASceneItem(A, h));
            elseif GAisa(A, 'bivector')
                if isempty(c)
                    c = 'g';
                end
                h = GAScene.drawhairydisk(PGAcast(A), c, origin(PGA));
                GAScene.additem(GASceneItem(A, h));
            elseif GAisa(A, 'trivector')
                if isempty(c)
                    c = 'b';
                end
                h = GAScene.drawhairyball(A, c, origin(PGA));
                GAScene.additem(GASceneItem(A, h));
            else
                error('Error is not a vector, bivector or trivector. OGA cannot draw it.');
            end

            GAScene.refreshdynamicitems();
        end

        function s = modelname(~)
            s = "OGA";
        end

        function R = cast(~, A)
            if isa(A, 'OGA')
                R = A;
            elseif isa(A, 'double')
                if GA.autoscalar
                    R = OGA(A);
                else 
                    error('Implicit conversion between a double and OGA is disabled. Run "help GA.autoscalar" for more information.')
                end
            else
                error(['Cannot implictly convert from ' class(A) ' to OGA'])
            end
        end

        function r = e1coeff(A)
            M = matrix(A);
            r = M(2);
        end

        function r = e2coeff(A)
            M = matrix(A);
            r = M(3);
        end

        function r = e3coeff(A)
            M = matrix(A);
            r = M(4);
        end

        function r = e12coeff(A)
            M = matrix(A);
            r = M(5);
        end

        function r = e13coeff(A)
            M = matrix(A);
            r = M(6);
        end

        function r = e23coeff(A)
            M = matrix(A);
            r = M(7);
        end

        function r = e123coeff(A)
            M = matrix(A);
            r = M(8);
        end
    end

end
classdef OldOGA
    % TODO: This is an old copy of PGA before the needed changes were made to upgrade dimensions.
    % Should probably just use this file for cross-referencing.
    properties
        m
    end

    % TODO: ensure return variables are all 'r' for consistency.
    % (perhaps not all 'r', but some naming convention.)
    
    % TODO: Add description of each method, including the "see also gable" stuff.

    %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%
    %           Settings           %
    %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%
    methods (Static = true)
        function val = innerProductType(newval)
            persistent currentval;
            
            % Default inner product is contraction
            if isempty(currentval)
                currentval = InnerProductType.Contraction;
            end

            if nargin >= 1
                if ismember(newval, enumeration('InnerProductType'))
                    currentval = newval;
                else
                    error('innerProductType must be an enumeration member of InnerProductType.')
                end
            end
            val = currentval;
        end

        function val = autoscalar(newval)
            persistent currentval;
            
            % Default inner product is contraction
            if isempty(currentval)
                currentval = true;
            end

            if nargin >= 1
                if islogical(newval)
                    currentval = newval;
                else
                    error('autoscalar must have value true or false.')
                end
            end
            val = currentval;
        end
    end
    

    % ZNOTE: TODO: consider merging this with the category below.

    methods (Access = private, Static = true)
        %ZNOTE: should actually be able to assume that returnPGA receives a PGA object. Thus, does not need to cast.
        % In fact, we should be able to make this non-static and private.
        % TODO: make the changes I noted above.
        function r = returnPGA(A)
            r = PGA(A);
            if PGAisa(r, 'scalar') && PGA.autoscalar
                r = r.m(1);
            end
        end
    end






































    
    %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%
    %    STATIC PRIVATE METHODS    %
    %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%

    % TODO: Make this section of methods private
    methods (Access = public, Static = true)

        % TODO: test this.
        function LineMesh(pts, c)
            if nargin == 1
                c = 'b';
            end
            M = size(pts, 1);
            N = size(pts, 2);
            for i=1:M
                for j=1:N
                    p = pts{i,j};
                    x(j) = p.m(2);
                    y(j) = p.m(3);
                    z(j) = p.m(4);
                end
                plot3(x,y,z,c);
            end
            for j=1:N
                for i=1:M
                    p = pts{i,j};
                    x(i) = p.m(2);
                    y(i) = p.m(3);
                    z(i) = p.m(4);
                end
                plot3(x,y,z,c);
            end
        end

        % TODO: test this.
        function a = gaarea(px,py)
            %gaarea(px,py): compute the area of a polygon.
            
            a = 0;
            for i = 1:length(px)-1
                a = a + px(i)*py(i+1) - px(i+1)*py(i);
            end
            a = a + px(length(px))*py(1) - px(1)*py(length(py));
            a = a/2;
        end

        function biarrow(iA, iB, c)
            % biarrow(A, B, c): draw an arrow from the head of B to the head of A in color c
            
            if nargin == 2
                c = 'b';
            end
            % ZNOTE: notice that we have to wrap PGAZ with PGA since PGA returns PGA/non-PGA based on user criteria
            %   If we were to instead make PGAZ private, we could not have this situation. Could be worth considering.
            A = PGA(PGAZ(iA));
            B = PGA(PGAZ(iB));
            
            % TODO: rephrase this in a better way
            if PGAisa(A, 'vector') == false || PGAisa(B,'vector') == false
               error('Can only draw arrow for vector');
            end
            
            % ZNOTE: need to convert this to PGA
            plot3([B.m(2) A.m(2)+B.m(2)], [B.m(3) A.m(3)+B.m(3)], [B.m(4) A.m(4)+B.m(4)], c);
            
            % Construct a vector perpendicular to A-B
            lA = sqrt(inner(A, A));
            % Numerical problems require us to extract the vector portion
            if abs(A.m(4)) < lA*.9
                %ZNOTE: This current fix is TEMPORARY!
                p1 = grade((PGA([0, 0, 0, 1, 0, 0, 0, 0])^A)*inverse(A), 1);
                %p1 = grade((e3^A)*inverse(A), 1);
            else
                %ZNOTE: This current fix is TEMPORARY!
                p1 = grade((PGA([0, 0, 1, 0, 0, 0, 0, 0])^A)*inverse(A), 1);
                %p1 = grade((e2^A)*inverse(A), 1);
            end
            hold on
            p2 = dual(p1^A);
            fprintf("First p2 check\n")
            display(p2)
            display(p1^A);
            display(dual(p1^A));
            lS = lA*.9;
            p1 = PGA((0.04*lS/sqrt(double(inner(p1, p1))))*p1);
            p2 = PGA((0.04*lS/sqrt(double(inner(p2, p2))))*p2);
            pA = A/lA*lS;
            % Cell array version
            t = (0:pi/4:2*pi);
            head = cell(2, length(t));
            % Pre add to reduce the cost below
            pAB = pA + B;
            AB = A + B;

            for i = 1:length(t)
                % ZNOTE: Used matrix here instead of .m, should this be the standard?
                % Work with matrices to reduce calls to GAExpand
                head{1, i} = PGA(sin(t(i))*matrix(p1) + cos(t(i))*matrix(p2) + matrix(pAB));
                head{2, i} = AB;
            end
            
            PGA.LineMesh(head, c);
        end

        function arrow(A, O, c)
            % arrowO(A,O,c): draw an arrow representing vector A in color c offset by O
            %  The c argument is optional; if not given, draw a blue arrow.

            A = PGA(A);
            O = PGA(O);

            % TODO: handle other nargin values.
            if nargin == 2
                c = 'b';
            end

            if ~isGrade(A, 1)
                error('Can only draw arrow for vector');
            end
            
            plot3([O.m(2) A.m(2)+O.m(2)], [O.m(3) A.m(3)+O.m(3)], [O.m(4) A.m(4)+O.m(4)], c);
            
            oldhold = ishold;
            hold on;
            
            % Construct a vector perpendicular to A
            lA = sqrt(inner(A, A));
            
            % Numerical problems require us to extra vector portion
            if abs(A.m(4)) < lA*.9
                %ZNOTE: Here temporarily using manual PGA(.) rather than e3, e2
                %p1 = grade((e3^A)*inverse(A), 1);
                p1 = grade((PGA([0, 0, 0, 1, 0, 0, 0, 0])^A)*inverse(A), 1);
            else
                %p1 = grade((e2^A)*inverse(A), 1);
                p1 = grade((PGA([0, 0, 1, 0, 0, 0, 0, 0])^A)*inverse(A), 1);
            end
            
            hold on
            p2 = dual(p1^A);
            lS = lA*.9;
            p1 = (0.04*lS/sqrt(double(inner(p1, p1))))*p1;
            p2 = (0.04*lS/sqrt(double(inner(p2, p2))))*p2;
            pA = A/lA*lS;

            % Cell array version
            t = (0:pi/4:2*pi);
            head = cell(2, length(t));

            % Pre-add to reduce the cost below
            pAO = pA + O;
            AO = A + O;

            for i = 1:length(t)
                % Work with matrices to avoid GAExpand calls
                head{1, i} = PGA(sin(t(i))*matrix(p1) + cos(t(i))*matrix(p2) + matrix(pAO));
                head{2, i} = AO;
            end
           
            PGA.LineMesh(head,c);
            if oldhold == 0
                hold off;
            end
        end

        function PGAMesh(pts, c)
            %PGAMesh(pts, c): Plot a mesh of GA vectors in color c (optional)

            if nargin == 1 
                c = .5;
            end
            M = size(pts, 1);
            N = size(pts, 2);
            
            for i = 1:M
                for j = 1:N
                p = pts{i,j};
                x(i, j) = p.m(2);
                y(i, j) = p.m(3);
                z(i, j) = p.m(4);
                C(i, j) = c;
                end
            end
            mesh(x, y, z, C);
        end

        function polygon(Poly, BV, O, c)
            % polygon(Poly,BV,O,c): draw a polygon
            %  Poly: a polygon struct
            %  BV: bivector plane in which to draw
            %  O: an offset from the origin
            %  c: the color of the polygon
            
            N = dual(BV);
            lA = abs(norm(N));
            
            if abs(N.m(4)) < lA*.9
               p1 = grade((e3^N)*inverse(N),1);
            else
               p1 = grade((e2^N)*inverse(N),1);
            end
            hold on
            p2 = dual(p1^N);
            
            p1 = (sqrt(lA/Poly.A)/sqrt(double(inner(p1,p1))))*p1;
            p2 = (sqrt(lA/Poly.A)/sqrt(double(inner(p2,p2))))*p2;
            
            % Cell array version
            for i=1:length(Poly.X)
                pts{i} = (Poly.X(i)-Poly.CX)*p1 + (Poly.Y(i)-Poly.CY)*p2 + O;
            end
            
            % Convert character color to RGB triple
            if isa(c,'char')
               if strncmp(c,'r',1)
                  c = [1 0 0];
               elseif strncmp(c,'g',1)
                  c = [0 1 0];
               elseif strncmp(c,'b',1)
                  c = [0 0 1];
               elseif strncmp(c,'c',1)
                  c = [0 1 1];
               elseif strncmp(c,'m',1)
                  c = [1 0 1];
               elseif strncmp(c,'y',1)
                  c = [1 1 0];
               elseif strncmp(c,'w',1)
                  c = [1 1 1];
               end
            end
            PGAPatch(pts,c);
        end

        %ZNOTE: Perhaps remove PGA prefix from this function and other static functions.
        function PGAPatch(pts, c)
            %PGAPatch(pts,c): Draw a patch in color c (optional)

            if nargin == 1
                c = 'y';
            end

            % Convert character color to RGB triple
            if isa(c, 'char')
                if strncmp(c, 'r', 1)
                    c = [1 0 0];
                elseif strncmp(c, 'g', 1)
                    c = [0 1 0];
                elseif strncmp(c, 'b', 1)
                    c = [0 0 1];
                elseif strncmp(c, 'c', 1)
                    c = [0 1 1];
                elseif strncmp(c, 'm', 1)
                    c = [1 0 1];
                elseif strncmp(c, 'y', 1)
                    c = [1 1 0];
                elseif strncmp(c, 'w', 1)
                    c = [1 1 1];
                end
            end
            for i = 1:length(pts)
                p = pts{i};
                x(i) = p.m(2);
                y(i) = p.m(3);
                z(i) = p.m(4);
            end
            patch(x, y, z, c);
        end

        function drawBivector(A, O, c)
            %drawBivector(A,O,c): Draw a bivector offset by vector O in color c.
            %  The c argument may be omitted.
            A = PGA(A);
            O = PGA(O);

            global GABVX;
            global GABVY;
            if ~strcmp(GAbvShape, 'default')
                PGA.arrow(0.2*dual(A), O, 'k');
            end
            if strcmp(GAbvShape, 'default')
                if nargin == 2
                    c = 'g';
                end
                % This line is a cheat - if your first object is 2D,
                %  then all subsequent draws will be 2D.  So draw
                %  a small 3D object first, hoping it won't be visible.
                PGA.arrow(0.001*dual(A),O,c);
                PGA.disk(A,O,c);
        
            elseif strcmp(GAbvShape, 'American')
                persistent gaAM;
                if isa(gaAM, 'double')
                    gaAM.X = [315 272 130 245 201 315 430 386 500 359 315];
                    gaAM.Y = [580 445 445 362 227 310 227 362 445 445 580];
                    gaAM.Z = [ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
                    gaAM.A = abs(gaarea(gaAM.X,gaAM.Y));
                    gaAM.CX = sum(gaAM.X)/length(gaAM.X);
                    gaAM.CY = sum(gaAM.Y)/length(gaAM.Y);
                end
                if nargin == 2
                    c = 'c';
                end
        
                polygon(gaAM, A, O, c);
            elseif strcmp(GAbvShape, 'Canadian')
                persistent gaCF;
                if isa(gaCF,'double')
                    gaCF.Y = [222 220 266 259 304 292 298 275 265 243 252 232 216    200   180   189   167   157   134   140   128   173   166   212   210];
                    gaCF.X = [ 3    43    38    51    90    96   126   121   135   110   161   157 188   157   161   110   135   121   126    96    90    51    38    43 3];
                    gaCF.Z = [ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
                    gaCF.A = abs(gaarea(gaCF.X,gaCF.Y));
                    gaCF.CX = sum(gaCF.X)/length(gaCF.X);
                    gaCF.CY = sum(gaCF.Y)/length(gaCF.Y);
                end
                if nargin == 2
                    c = 'r';
                end
        
                polygon(gaCF, A, O, c);
            elseif strcmp(GAbvShape, 'Dutch')
                persistent gaWM;
                if isa(gaWM,'double')
                    gaWM.Y = [134 161 137 121 145 163 171 125 124 171 165 143 119 91 72 41 22 16 59 58 19 26 48 48 18 64];
                    gaWM.X = -1*[229 201 179 119 138 128 108  86  74  20   5   0  38 17 30 5  16 34 57 82 120 136 144 177 209 234];
                    gaWM.Z = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
            
                    gaWM.A = abs(gaarea(gaWM.X,gaWM.Y));
                    gaWM.CX = sum(gaWM.X)/length(gaWM.X);
                    gaWM.CY = sum(gaWM.Y)/length(gaWM.Y);
                end
                if nargin == 2
                    c = 'b';
                end
                polygon(gaWM, A, O, c);
            elseif strcmp(GAbvShape, 'Clifford')
                persistent gaCL;
                if isa(gaCL,'double')
                    gaCL.Y = [ 5850 5841 5700 5325 5100 4875 5100 5100 4950 5100 5550 5400 5550 5850 5550 5400 5550 5100 5100 5325 5700 5850 5925 5850 5710 5448 5280 5448 5719 5700 5775 5700 5775 5550 5250 5250 5700 5250 5175 5250 5175 5175 5250 5700 5250 5100 4950 5100 4800 4725 4650 4725 4425 4500 4500 4650 4950 4950 4725 4725 4875 4800 4800 5025 5025 5025 5175 5550 5175 5250 5550 5250 4875 4875 5325 5400 5550 5550 5550 5850 5550 5625 5550 5550 5700 5700 5850 5925 5775 5850 5850];
                    gaCL.X = -1*[ 4275 4226 4275 4350 4275 4275 4275 4125 4050 4125 3975 3675 3375 3300 3375 3675 3975 4125 4275 4350 4275 4200 3975 3300 3020 3235 3235 3226 3001 2775 2550 2475 2475 2250 2325 2475 2475 2475 2550 2700 3075 3225 3225 2850 3225 3525 3900 3525 3600 3600 3825 3600 3675 3675 3825 3825 3900 4050 4050 4200 4275 4950 5700 5700 4875 6000 6000 6000 6000 6450 6450 6450 6600 6675 6600 6675 6600 6450 6000 5775 6000 6225 6375 6450 6450 6525 6525 6300 4950 4875 4275];
                    gaCL.A = abs(gaarea(gaCL.X,gaCL.Y));
                    gaCL.CX = sum(gaCL.X)/length(gaCL.X);
                    gaCL.CY = sum(gaCL.Y)/length(gaCL.Y);
                end
                if nargin == 2
                    c = 'b';
                end
                polygon(gaCL, A, O, c);
                
            elseif strcmp(GAbvShape, 'CLIFFORD')
                % Make Clifford Red
                c = 'r';
                persistent gaCBRD;
                if isa(gaCBRD,'double')
                    gaCBRD.Y = [92    83    77    70    65    55    51    53    66    66    64    41 31    10    10    16    29    25    31    27    26    29    35    32 35    33    35    31    35    33    35    48    51    50    41    44 48    48    44    48    48    44    58    58    61    66    68    64 65    68    77    80    86    87    80    74    77    77    79    89 93    88    86    82    79    76    79    82    86    88    93    93 92];
                    gaCBRD.X = -1*[35    35    30    29    35    32    36    40    40    39    45    43 45    21    30    41    48    56    66    72    84    86    83    80 73    70    67    66    67    70    73    80    78    75    68    66 60    52    48    52    60    66    69    84    86    86    84    81 68    66    67    71    70    67    62    62    56    49    46    48 43    43    45    43    44    40    44    43    45    43    43    39 35];
                    gaCBRD.Z = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
                    %gaCBRD.A = abs(gaarea(gaCBRD.X,gaCBRD.X));
                    % Make Clifford Big
                    gaCBRD.A = 100;
                    gaCBRD.CX = sum(gaCBRD.X)/length(gaCBRD.X);
                    gaCBRD.CY = sum(gaCBRD.Y)/length(gaCBRD.Y);
                end
                polygon(gaCBRD, A, O, c);
            elseif strcmp(GAbvShape, 'UserDefined')
                if nargin == 2
                    c = 'c';
                end
                gaUP.X = GABVX;
                gaUP.Y = GABVY;
                gaUP.Z = zeros(size(GABVX));
                gaUP.A = abs(gaarea(gaUP.X,gaUP.Y));
                gaUP.CX = sum(gaUP.X)/length(gaUP.X);
                gaUP.CY = sum(gaUP.Y)/length(gaUP.Y);
                polygon(gaUP, A, O, c);
            else
                error('Unknown GAbvShape');
            end
        end

        % ZNOTE: It's not clear what O is, conceptually. Perhaps a different data structure would be better suited
        % Also, extracting the vector terms using .(2), .(3), .(4) is opaque. Perhaps we should consider writing a function
        % to extract the vector/bivector etc. portions of multivectors more explictly (this could help with upgrading the code
        % as well)
        function drawTrivector(tv, O, c)
            %drawTrivector(tv, O, c): draw a trivector
            tv = PGA(tv);
            O = PGA(O);

            if nargin == 2
                c = 'r';
            end
            v = abs(tv.m(8));
            r = (v*3/4/pi)^(1./3.);
        
            [X,Y,Z]=sphere(8);
            X = r*X;
            Y = r*Y;
            Z = r*Z;
            plot3(X+O.m(2),Y+O.m(3),Z+O.m(4),c);
            hold on
            for i=1:size(X,1)
                for j=1:size(X,2)
                    x(j) = X(i,j);
                    y(j) = Y(i,j);
                    z(j) = Z(i,j);
                end
                plot3(x+O.m(2),y+O.m(3),z+O.m(4),c);
            end
            if tv.m(8) > 0
                d = 1.2;
            else
                d = .85;
            end

            for i=1:size(X,1)
                for j=mod(i,2)+1:2:size(X,2)
                    xh(1,j) = X(i,j);
                    xh(2,j) = d*X(i,j);
                    yh(1,j) = Y(i,j);
                    yh(2,j) = d*Y(i,j);
                    zh(1,j) = Z(i,j);
                    zh(2,j) = d*Z(i,j);
                end
                plot3(xh+O.m(2),yh+O.m(3),zh+O.m(4),c);
            end
        end

        function r = DrawProduct(iA, iB, p)
            %DrawProduct(A,B,p): draw a product.
            % A,B: multivectors from which to form the product.
            % p: string declaring type of product.
            
            clf;
            A = PGA(PGAZ(iA));
            B = PGA(PGAZ(iB));
            
            if strcmp(p,'inner')
                r = inner(A,B);
            elseif strcmp(p, 'wedge')
                r = A^B;
            elseif strcmp(p, 'geometric')
                r = A*B;
            elseif strcmp(p, 'spinor')
                r = inverse(A)*B*A;
            end
            
            set(gcf,'Name',p);
            
            if PGAisa(A, 'multivector') || PGAisa(B, 'multivector')
                subplot(1,3,1);
                draw(A);
                if isa(A, 'double') || PGAisa(A, 'double')
                    a = [0 0 0 0 0 0];
                    axis('equal')
                else
                    a = axis;
                end
                subplot(1,3,2);
                draw(B);
                if isa(B, 'double') || PGAisa(B, 'double')
                    b = [0 0 0 0 0 0];
                    axis('equal')
                else
                  b = axis;
                end
                subplot(1,3,3);
                draw(r);
                if isa(r, 'double') || PGAisa(r, 'double')
                    c = [0 0 0 0 0 0];
                    axis('equal')
                else
                    c = axis;
                end
                for i=1:2:5
                    d(i) = min([a(i),b(i),c(i)]);
                    d(i+1) = max([a(i+1),b(i+1),c(i+1)]);
                end
                subplot(1,3,1);axis(d);
                subplot(1,3,2);axis(d);
                subplot(1,3,3);axis(d);
            else
                subplot(1,2,1);
                draw(A,'b');
                draw(B,'g');
                if isa(A,'double') && PGAisa(B,'double')
                    a = [0 0 0 0 0 0];
                    axis('equal')
                else
                    a = axis;
                end
                subplot(1,2,2);
                draw(r, 'r');
                if isa(r,'double') || PGAisa(r,'double')
                    b = [0 0 0 0 0 0];
                    axis('equal')
                else
                    b = axis;
                end
                for i = 1:2:5
                    d(i) = min([a(i),b(i)]);
                    d(i+1) = max([a(i+1),b(i+1)]);
                end
                subplot(1,2,1);
                axis(d);
                subplot(1,2,2);
                axis(d);
            end
        end 
        
        function disk(A,O,c)
            %disk: draw a disk to represent bivector A offset by O in color c.
            A = PGA(A);
            O = PGA(O);

            if ~isGrade(A, 2)
                error('Can only draw disk for bivector');
            end
        
            dA = dual(A);
    
            lA = abs(norm(dA));
            if abs(dA.m(4)) < lA*.9
                % ZNOTE: TODO: here using PGA(.) instead of e3, e2
                % p1 = grade((e3^dA)*inverse(dA),1);
                p1 = grade((PGA([0, 0, 0, 1, 0, 0, 0, 0])^dA)*inverse(dA),1);
            else
                %p1 = grade((e2^dA)*inverse(dA),1);
                p1 = grade((PGA([0, 0, 1, 0, 0, 0, 0, 0])^dA)*inverse(dA),1);
            end

            hold on
            p2 = dual(p1^dA);
    
            p1 = (sqrt(lA/pi)/sqrt(abs(double(inner(p1,p1)))))*p1;
            p2 = (sqrt(lA/pi)/sqrt(abs(double(inner(p2,p2)))))*p2;
            % Cell array version
            t = (0:pi/12:2*pi);
            for i=1:length(t)
                % Speed it up a little by using matrices directly
                pts{i} = PGA(sin(t(i))*matrix(p1) + cos(t(i))*matrix(p2) + matrix(O));
            end
    
            PGA.PGAPatch(pts,c);
            for i=1:4:length(t)-1
                p{1} = pts{i};
                p{2} = 1.1*(pts{i}-O) + dual((pts{i}-O)^dA)/abs(norm(dA))/4+O;
                PGADrawPolyline({p{2},p{1}},'k');
            end
        end
    end



















































    %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%
    % PRIVATE (NON-STATIC) METHODS %
    %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%
    
    methods (Access = private)

        function r = PGAisa_(A, t)
            nzm = A.m ~= 0;
            if strcmp(t, 'double') || strcmp(t, 'scalar') 
                r = sum(nzm(2:8)) == 0;
            elseif strcmp(t,'vector') 
                r = sum(nzm(5:8)) == 0 & nzm(1) == 0;
            elseif strcmp(t,'bivector') 
                r = sum(nzm(1:4)) == 0 & nzm(8) == 0;
            elseif strcmp(t,'trivector') 
                r = sum(nzm(1:7)) == 0;
            elseif strcmp(t,'multivector')
                r = sum( [sum(nzm(1)) sum(nzm(2:4)) sum(nzm(5:7)) nzm(8)] ~= 0) > 1;
            else
                r = false;
            end 
        end

        % ZNOTE: Having a private counterpart to this method is not necessary. Do it anyway?
        function r = double_(A)
            if PGAisa_(A, 'scalar')
                r = A.m(1);
            else
                error('Can only convert a scalar PGA object to a double.');
            end
        end

        % ***** Functions for adding and subtracting PGA objects *****

        function r = plus_(A, B)
            r = PGA(A.m + B.m);
        end

        function r = minus_(A, B)
            r = PGA(A.m - B.m);
        end


        % ***** Geometric Product Stuff *****

        % TODO: rename to PGAsignature
        % TODO: Consider removing this functionality, perhaps?
        function G = GAsignature(Sign1,Sign2,Sign3)
            persistent GAsignature3;
            persistent GAsignature2;
            persistent GAsignature1;
            if nargin == 3
               GAsignature3=Sign3;
               GAsignature2=Sign2;
               GAsignature1=Sign1;
            end
            G = [prod(GAsignature1) ,prod(GAsignature2), prod(GAsignature3)];
        end

        % TODO: lower case E? or not?
        % TODO: make private
        function r = PGAExpand(A)
            S = GAsignature;
            r = [A.m(1)  S(1)*A.m(2)  S(2)*A.m(3)  S(3)*A.m(4) -S(1)*S(2)*A.m(5) -S(2)*S(3)*A.m(6) -S(1)*S(3)*A.m(7) -S(1)*S(2)*S(3)*A.m(8);
                A.m(2)     A.m(1)    S(2)*A.m(5) -S(3)*A.m(7)   -S(2)*A.m(3)    -S(2)*S(3)*A.m(8)    S(3)*A.m(4)     -S(2)*S(3)*A.m(6);
                A.m(3) -S(1)*A.m(5)     A.m(1)    S(3)*A.m(6)    S(1)*A.m(2)      -S(3)*A.m(4)    -S(3)*S(1)*A.m(8)  -S(1)*S(3)*A.m(7);
                A.m(4)  S(1)*A.m(7) -S(2)*A.m(6)     A.m(1)   -S(1)*S(2)*A.m(8)    S(2)*A.m(3)      -S(1)*A.m(2)     -S(1)*S(2)*A.m(5);
                A.m(5)    -A.m(3)       A.m(2)    S(3)*A.m(8)       A.m(1)         S(3)*A.m(7)      -S(3)*A.m(6)         S(3)*A.m(4);
                A.m(6)  S(1)*A.m(8)    -A.m(4)       A.m(3)     -S(1)*A.m(7)          A.m(1)         S(1)*A.m(5)         S(1)*A.m(2);
                A.m(7)     A.m(4)    S(2)*A.m(8)    -A.m(2)      S(2)*A.m(6)      -S(2)*A.m(5)          A.m(1)           S(2)*A.m(3);
                A.m(8)     A.m(6)       A.m(7)       A.m(5)         A.m(4)            A.m(2)            A.m(3)              A.m(1)];
        end

        function r = product_(A, B)
            A_expanded = PGAExpand(A);
            r = PGA(A_expanded*B.m);
        end

        % ***** The Outer Product *****

        function r = outer_(A, B)
                M = [A.m(1)    0       0       0       0       0       0       0;
                     A.m(2)  A.m(1)    0       0       0       0       0       0;
                     A.m(3)    0     A.m(1)    0       0       0       0       0;
                     A.m(4)    0       0     A.m(1)    0       0       0       0;
                     A.m(5) -A.m(3)  A.m(2)    0     A.m(1)    0       0       0;
                     A.m(6)    0    -A.m(4)  A.m(3)    0     A.m(1)    0       0;
                     A.m(7)  A.m(4)    0    -A.m(2)    0       0     A.m(1)    0;
                     A.m(8)  A.m(6)  A.m(7)  A.m(5)  A.m(4)  A.m(2)  A.m(3)  A.m(1)];
                
                r = PGA(M*B.m);
        end

        % ***** The Inner Product *****
        function r = inner_(A, B)
            switch(PGA.innerProductType)
                case InnerProductType.Contraction
                    r = contraction_(A,B);
                case InnerProductType.Hestenes
                    r = innerH_(A,B);
                case InnerProductType.HestenesModifiedForScalars
                    r = innerS_(A,B);
                otherwise
                    error('Unknown inner type');
            end
        end

        function r = contraction_(A, B)
            S=GAsignature;
            M = [A.m(1) S(1)*A.m(2) S(2)*A.m(3) S(3)*A.m(4) -S(1)*S(2)*A.m(5) -S(2)*S(3)*A.m(6) -S(3)*S(1)*A.m(7) -S(1)*S(2)*S(3)*A.m(8);
	            0      A.m(1)        0             0       -S(2)*A.m(3)            0            S(3)*A.m(4)    -S(2)*S(3)*A.m(6);
	            0        0         A.m(1)          0        S(1)*A.m(2)       -S(3)*A.m(4)           0         -S(1)*S(3)*A.m(7);
	            0        0           0           A.m(1)          0             S(2)*A.m(3)     -S(1)*A.m(2)    -S(1)*S(2)*A.m(5);
	            0        0           0             0           A.m(1)              0                 0           S(3)*A.m(4);
	            0        0           0             0             0               A.m(1)              0           S(1)*A.m(2);
	            0        0           0             0             0                 0               A.m(1)        S(2)*A.m(3);
	            0        0           0             0             0                 0                 0             A.m(1)];
            r = PGA(M*B.m);
        end

        function r = innerH_(A, B)
            S=GAsignature;  
            Mga = [ 0  S(1)*A.m(2)  S(2)*A.m(3)  S(3)*A.m(4) -S(1)*S(2)*A.m(5) -S(2)*S(3)*A.m(6) -S(1)*S(3)*A.m(7) -S(1)*S(2)*S(3)*A.m(8);
	                0       0       S(2)*A.m(5) -S(3)*A.m(7)   -S(2)*A.m(3)    -S(2)*S(3)*A.m(8)    S(3)*A.m(4)     -S(2)*S(3)*A.m(6);
	                0 -S(1)*A.m(5)       0       S(3)*A.m(6)    S(1)*A.m(2)      -S(3)*A.m(4)    -S(3)*S(1)*A.m(8)  -S(1)*S(3)*A.m(7);
	                0  S(1)*A.m(7) -S(2)*A.m(6)       0      -S(1)*S(2)*A.m(8)    S(2)*A.m(3)      -S(1)*A.m(2)     -S(1)*S(2)*A.m(5);
	                0       0            0       S(3)*A.m(8)         0                 0                 0              S(3)*A.m(4);
	                0  S(1)*A.m(8)       0            0              0                 0                 0              S(1)*A.m(2);
	                0       0       S(2)*A.m(8)       0              0                 0                 0              S(2)*A.m(3);
	                0       0            0            0              0                 0                 0                   0];
            r = PGA(Mga*B.m);
        end

        function r = innerS_(A, B)
            S = GAsignature;  
            Mga = [ A.m(1)  S(1)*A.m(2)  S(2)*A.m(3)  S(3)*A.m(4) -S(1)*S(2)*A.m(5) -S(2)*S(3)*A.m(6) -S(1)*S(3)*A.m(7) -S(1)*S(2)*S(3)*A.m(8);
                A.m(2)       A.m(1)  S(2)*A.m(5) -S(3)*A.m(7)   -S(2)*A.m(3)    -S(2)*S(3)*A.m(8)    S(3)*A.m(4)     -S(2)*S(3)*A.m(6);
                A.m(3) -S(1)*A.m(5)       A.m(1)  S(3)*A.m(6)    S(1)*A.m(2)      -S(3)*A.m(4)    -S(3)*S(1)*A.m(8)  -S(1)*S(3)*A.m(7);
                A.m(4)  S(1)*A.m(7) -S(2)*A.m(6)       A.m(1) -S(1)*S(2)*A.m(8)    S(2)*A.m(3)      -S(1)*A.m(2)     -S(1)*S(2)*A.m(5);
                A.m(5)       0            0       S(3)*A.m(8)         A.m(1)            0                 0              S(3)*A.m(4);
                A.m(6)  S(1)*A.m(8)       0            0              0                 A.m(1)            0              S(1)*A.m(2);
                A.m(7)       0       S(2)*A.m(8)       0              0                 0                 A.m(1)         S(2)*A.m(3);
                A.m(8)       0            0            0              0                 0                 0                   A.m(1)];
            r = PGA(Mga*B.m);
        end

        % ***** Inverse and Division *****

        function r = smallinverse_(A)
            S = prod(GAsignature);
            P = A.m(1)*A.m(1) + A.m(8)*A.m(8)*S;
            if P==0
               disp('smallinverse: This multivector does not have an inverse.');
               r = PGA(0);
            else
               r = PGA([A.m(1)/P; 0; 0; 0; 0; 0; 0; -A.m(8)/P]);
            end
        end

        function r = inverse_(A)
            % G is the Clifford conjugate of A
            G = PGA([A.m(1); -A.m(2); -A.m(3); -A.m(4); - A.m(5); -A.m(6);  -A.m(7);  A.m(8)]);
            V = product_(A, G);	% V = AG
            % If AG is scalar, then inverse(AG)=1/(AG) and inverse is cheap
            %  V is of the form Scalar+Pseudoscalar, so to test if it is
            %  a scalar, we just check the pseudoscalar coefficient.
            if V.m(8) == 0
                Den = V.m(1);
                if Den == 0
                    disp('This multivector does not have an inverse.');
                    r = PGA(0);
                else 
                    r = PGA(1/Den).product_(G);
                end
            else
                % Conceptually we want
                %        Mi=PGAproduct(smallinverse(V),G);
                % But the following is more efficient
                Vi = smallinverse_(V);
                S = prod(GAsignature);
                Mi = G.product_(PGA(Vi.m(1))) - PGA(S*Vi.m(8)).product_(PGA(dual_(G)));

                r = PGA(Mi);
            end
        end     

        % ***** Display *****

        function s = charF(m, f)
            L=f.m;
            R=m.m;
            p=L*R;
            n1=f.b1;
            n2=f.b2;
            n3=f.b3;
            pl = ' ';
            s = '    ';
            if p(1) ~= 0
                s = [s pl num2str(p(1))];
                pl = ' + ';
            end
            if p(2) ~= 0
                if p(2) == 1
                    s = [s pl n1];
                else
                    s = [s pl num2str(p(2)) '*' n1];
                end
                pl = ' + ';
            end
            if p(3) ~= 0
                if p(3) == 1
                    s = [s pl n2];
                else
                    s = [s pl num2str(p(3)) '*' n2];
                end
                pl = ' + ';
            end
            if p(4) ~= 0
            if p(4) == 1
                s = [s pl n3];
            else
                s = [s pl num2str(p(4)) '*' n3];
            end
            pl = ' + ';
            end

            if p(5) ~= 0
            if p(5) == 1
                s = [s pl n1 '^' n2];
            else
                s = [s pl num2str(p(5)) '*' n1 '^' n2];
            end
            pl = ' + ';
            end
            if p(6) ~= 0
                if p(6) == 1
                    s = [s pl n2 '^' n3];
                else
                    s = [s pl num2str(p(6)) '*' n2 '^' n3];
                end
                pl = ' + ';
            end
            if p(7) ~= 0
                if p(7) == 1
                    s = [s pl n3 '^' n1];
                else
                    s = [s pl num2str(p(7)) '*' n3 '^' n1];
                end
                pl = ' + ';
            end

            if p(8) ~= 0
                if p(8) == 1
                    s = [s pl n1 '^' n2 '^' n3];
                else
                    s = [s pl num2str(p(8)) '*' n1 '^' n2 '^' n3];
                end
                pl = ' + ';
            end

            if strcmp(pl, ' ')
            s = '     0';
            end
        end

        function s = char(p)
            if isa(OFrame,'struct')
                s = charF(p,OFrame);
            else
                pl = ' ';
                s = '    ';
                if p.m(1) ~= 0
                    s = [s pl num2str(p.m(1))];
                    pl = ' + ';
                end
                if p.m(2) ~= 0
                    if p.m(2) == 1
                        s = [s pl 'e1'];
                    else
                        s = [s pl num2str(p.m(2)) '*e1'];
                    end
                    pl = ' + ';
                end
                if p.m(3) ~= 0
                    if p.m(3) == 1
                        s = [s pl 'e2'];
                    else
                        s = [s pl num2str(p.m(3)) '*e2'];
                    end
                    pl = ' + ';
                end
                if p.m(4) ~= 0
                    if p.m(4) == 1
                        s = [s pl 'e3'];
                    else
                        s = [s pl num2str(p.m(4)) '*e3'];
                    end
                    pl = ' + ';
                end
                
                if p.m(5) ~= 0
                    if p.m(5) == 1
                        s = [s pl 'e1^e2'];
                    else
                        s = [s pl num2str(p.m(5)) '*e1^e2'];
                    end
                    pl = ' + ';
                end
                if p.m(6) ~= 0
                    if p.m(6) == 1
                        s = [s pl 'e2^e3'];
                    else
                        s = [s pl num2str(p.m(6)) '*e2^e3'];
                    end
                    pl = ' + ';
                end
                if p.m(7) ~= 0
                    if p.m(7) == 1
                        s = [s pl 'e3^e1'];
                    else
                        s = [s pl num2str(p.m(7)) '*e3^e1'];
                    end
                    pl = ' + ';
                end
                
                if p.m(8) ~= 0
                    if p.m(8) == 1
                        s = [s pl 'I3'];
                    else
                        s = [s pl num2str(p.m(8)) '*I3'];
                    end
                    pl = ' + ';
                end
                
                if strcmp(pl, ' ')
                  s = '     0';
                end
            end
        end

        % TODO: probably should rename this function.
        % TODO: perhaps depricate this?
        function s = newchar(p)
            %newchar(p): convert a GA object to a string proclaiming its type.

            pl = ' ';
            s = '    ';
            if PGAisa(p,'scalar')
              s = [s 'scalar'];
            elseif PGAisa(p, 'vector')
              s = [s 'vector'];
            elseif PGAisa(p, 'bivector')
              s = [s 'bivector'];
            elseif PGAisa(p, 'trivector')
              s = [s 'trivector'];
            else
              s = [s 'multivector'];
            end
        end

        % ***** Norms *****

        function r = norm_(A)
            B = bilinear_(A, A);
            if B > 0
                r = sqrt(B);
            else
                r = sqrt(-B);
            end
        end

        % ***** Orthogonalize and Frame *****

        % TODO: test this
        function N = Orthogonalize(normal,basis,CLsize)
            %Orthogonalize(n,b,s): Gives a basis orthogonal to the normal.
            %First argument is the normal vector, second is the basis matrix
            %third is the size of a PGA (optional).
            if nargin < 3
                CLsize = 8;
            end
            N = [];  
            ni = normal;   
            Imax = size(basis,2);
            for i = Imax:-1:1
                Ii = inverse(ni);
                for j = 1:i
                    P = PGAproduct(PGAouter(PGA(basis(1:CLsize,j)),ni),Ii);
                    R(1:CLsize,j) = P.m;
                end
                T = R(1:CLsize,1);
                ni = PGA(T)
                v = norm(ni);
                if v == 0
                    M = T;
                else
                    M = (1/v)*T;
                end
                N = [N,M];
                basis = basis(1:CLsize, 2:size(basis, 2));
            end
        end

        % TODO: test this
        function F = Frame(fname, v1, v2, v3, n1, n2, n3)
            %Frame(fname,v1,v2,v3,n1,n2,n3): Creates a frame for PGA basis vectors 
            %Takes a name, and three linearly independent vectors (of type PGA) 
            % and optionally the three vectors can be named, 
            %If they are not named, then the names are the frame name followed by 
            % 1,2, or 3 respectively.
            if nargin < 4
                disp('The frame needs a name, and three linearly independent vectors.');
                F.m = eye(8);
                F.b1 = 'e1';
                F.b2 = 'e2';
                F.b3 = 'e3';
            else
                if nargin < 7
                    n1 = [fname, '1'];
                    n2 = [fname, '2'];
                    n3 = [fname, '3'];
                else
                    %TODO: why does this else exist?
                end

                %TODO: Change this to casting to PGA
                if isa(v1,'PGA') && isa(v2,'PGA') && isa(v3,'PGA')
                    E = v1^v2^v3;
                    EM = E.m;
                    W = EM(8);
                    if W == 0
                        disp('The frame needs three linearly independent vectors.');
                        F.m=eye(8);
                        F.b1='e1';
                        F.b2='e2';
                        F.b3='e3';
                    else
                        % ZNOTE: TODO: This used to be Edual. Probably necessary for Frame.
                        r1 = (1/W)*dual(v2^v3);
                        r2 = (1/W)*dual(v3^v1);
                        r3 = (1/W)*dual(v1^v2);
                        %Bivectors
                        r12 = r1^r2;
                        r23 = r2^r3;
                        r31 = r3^r1;
                        G = [1,0,0,0,0,0,0,0; 
                            r1.m'; 
                            r2.m'; 
                            r3.m';
                            r12.m';
                            r23.m';
                            r31.m';
                            0,0,0,0,0,0,0,1/W];
                        F.m = G;
                        F.b1 = n1;
                        F.b2 = n2;
                        F.b3 = n3;
                    end
                else
                    disp('The frame needs three linearly independent vectors.');
                    F.m = eye(8);
                    F.b1 = 'e1';
                    F.b2 = 'e2';
                    F.b3 = 'e3';
                end
            end
        end



        % ***** Equality and Inequality *****

        function r = eq_(A, B)
            % TODO: Perhaps make tolerance adjustable? Or also make it so eeq can be made default?
            tol = 1e-15;

            r = norm(A - B) < tol;
        end

        function r = eeq_(A, B)
            r = all(A.m == B.m);
        end

        function r = ne_(A, B)
            r = ~eq(A, B);
        end
        
        % ***** Dual and Reverse*****

        function r = dual_(A)
            %The dual is calculated by right multiplication by the inverse of the
            % pseudoscalar. In this case we get a negative sign from the grade of 
            % the pseudoscalar and the other term from the signature.
            S = GAsignature;
            % Conceptually we want
            %    T=GAproduct(A,GA([0;0;0;0;0;0;0;-prod(S)]));
            % but the code is more efficient as
            r = PGA([A.m(8); S(1)*A.m(6); S(2)*A.m(7); S(3)*A.m(5); -S(1)*S(2)*A.m(4); -S(2)*S(3)*A.m(2); -S(3)*S(1)*A.m(3); -S(1)*S(2)*S(3)*A.m(1)]);
        end

        function r = reverse_(A)
            r = PGA([A.m(1);A.m(2);A.m(3);A.m(4);-A.m(5);-A.m(6);-A.m(7);-A.m(8)]);
        end

        % ***** Meet and Join *****

        function r = PGAZ(A)
            % PGAZ(A): Set to zeros any element of the multivector A less than 1e-15
            r = PGA(A);
            for i=1:8
              if abs(r.m(i)) < 1e-15
                 r.m(i) = 0;
              end
            end
        end

        function r = join_(A, B)
            % Approach: We will exploit the structure of our 3D space.
            %  In particular, the wedge product will often (but not) always
            %  compute what we want.  When it fails, the dual product may
            %  compute what we want.  Together, the two of them only fail
            %  when both objects are vectors or bivectors, and the join
            %  is degenerate.  In such cases, the join will either be the
            %  bivector, or if both arguments are vectors, then its the
            %  vector.
            % Note: this approach works because except in degenerate cases,
            %  if a and b have a common subspace, then dual(a) and dual(b) 
            %  do not have a common subspace.  Thus, we can not easily
            %  generalize this approach to higher dimensional spaces.

            a = PGAZ(A);
            b = PGAZ(B);
            if PGAisa(a,'multivector') || PGAisa(b,'multivector')
                error('join: both arguments must be blades.');
            end
            p = PGAZ(a^b);
            if p ~= 0
                r = p;
            else
                R = PGAZ(contraction(dual(b),a));
                if R ~= 0
                    r = norm(R)*(a/R)^b;
                else
                    if PGAisa(a,'bivector')
                        r = norm(b)*a;	
                    else
                        r = norm(a)*b;	
                    end
                end
            end

            r = PGA(r);
        end

        function r = meet_(A, B)
            a = PGAZ(A);
            b = PGAZ(B);

            if PGAisa(a,'multivector') || PGAisa(b,'multivector') 
                error('meet: all arguments must be blades.');
            end

            r = contraction_(b/join(a,b), a);
        end

        % ***** sLog, conjugate, gradeInvolution *****

        function r = sLog_(spinor)
            % Takes logarithm of spinor
            % TODO: Spinor verification
            
            plane = unit(grade_(spinor, 2));
            r = log(norm(spinor)) + plane*atan2(double(grade(spinor/plane, 0)), double(grade(spinor, 0)));
        end

        function r = conjugate_(A)
            %conjugate(A): compute the Clifford conjugate of a multivector.
            
            r = PGA([A.m(1); -A.m(2); -A.m(3); -A.m(4); - A.m(5); -A.m(6);  -A.m(7);  A.m(8)]);
        end

        function r = gradeinvolution_(A)
            %gradeinvolution(A): returns the grade involution of a multivector.
            r = PGA([A.m(1);-A.m(2);-A.m(3);-A.m(4);A.m(5);A.m(6);A.m(7);-A.m(8)]);
        end

        % ***** Grade, blade *****

        function r = grade_(A, n)
            % grade(A, n): return the part of an object of a particular grade.
            %  Return the part of M of grade n.
            %  If n is omitted, return the grade of M (-1 if M is of mixed grade) 
            
            if nargin == 1
                if A.m(1) ~= 0
                    if sum(abs(A.m(2:8))) == 0
                        r = 0;
                    else
                        r = -1;
                    end
                elseif sum(abs(A.m(2:4))) ~= 0
                    if sum(abs(A.m(5:8))) == 0
                        r = 1;
                    else
                        r = -1;
                    end
                elseif sum(abs(A.m(5:7))) ~= 0
                    if A.m(8) == 0
                        r = 2;
                    else
                        r = -1;
                    end
                elseif A.m(8) ~= 0
                    r = 3;
                else
                    r = -1;
                end
            else
                if n == 0
                    r = PGA([A.m(1);0;0;0;0;0;0;0]);
                elseif n == 1
                    r = PGA([0;A.m(2);A.m(3);A.m(4);0;0;0;0]);
                elseif n == 2
                    r = PGA([0;0;0;0;A.m(5);A.m(6);A.m(7);0]);
                elseif n == 3
                    r = PGA([0;0;0;0;0;0;0;A.m(8)]);
                else
                    r = PGA(0);
                end
            end
        end

        function r = isGrade_(A, g)
            r = false;

            if g==0
                if sum(abs(A.m(2:8))) == 0
                    r = true;
                end
            elseif g==1
                if sum(abs([A.m(1);A.m(5:8)])) == 0
                    r = true;
                end
            elseif g==2
                if sum(abs([A.m(1:4);A.m(8)])) == 0
                    r = true;
                end
            elseif g==3
                if sum(abs(A.m(1:7))) == 0
                    r = true;
                end
            elseif g==-1
                z = A.m == 0;
                z0 = z(1);
                z1 = sum(z(2:4))~=0;
                z2 = sum(z(5:7))~=0;
                z3 = z(8);
                % Note that the test treats 0 is a multivector!
                if z0+z1+z2+z3 ~= 1
                    r = true;
                end
            else
                error('isGrade: invalid grade.');
            end
        end
            

        function r = unit(A)
            if isGrade(A,0)
                r = unit(A.m(1));
            elseif isGrade(A,1)
                s = sqrt(A.m(2)*A.m(2)+A.m(3)*A.m(3)+A.m(4)*A.m(4));
                r = A/s;
            elseif isGrade(A,2)
                s = sqrt(A.m(5)*A.m(5)+A.m(6)*A.m(6)+A.m(7)*A.m(7));
                r = A/s;
            elseif isGrade(A,3)
                r = A/abs(A.m(8));
            else
                error('Unit can only be applied to blades.');
            end
        end

        function r = blade(A)
            % blade(A) : return a blade made from the largest portion of a multivector.
            A = PGA(A);

            s(1) = abs(A.m(1));
            s(2) = sqrt(sum(abs(A.m(2:4))));
            s(3) = sqrt(sum(abs(A.m(5:7))));
            s(4) = abs(A.m(8));
            if s(1)>s(2) && s(1)>s(3) && s(1)>s(4)
              r = PGA.returnPGA(A.m(1));
            elseif s(2)>s(3) && s(2)>s(4)
              r = PGA.returnPGA([0; A.m(2); A.m(3); A.m(4); 0; 0; 0; 0]);
            elseif s(3)>s(4)
              r = PGA.returnPGA([0; 0; 0; 0; A.m(5); A.m(6); A.m(7); 0]);
            else
              r = PGA.returnPGA([0; 0; 0; 0; 0; 0; 0; A.m(8)]);
            end
        end

        function r = gexp_(A)
            %gexp(A): Computes the geometric product exponential of a multivector.
            M = PGAExpand(A);
            E = expm(M);
            r = PGA(E(1:8,1));
        end

        function r = wexp_(A)
            %GAwexp(A): Gives the wedge product exponential of a GA objects.
        
            M =[A.m(1)    0       0       0       0       0       0       0;
            A.m(2)  A.m(1)    0       0       0       0       0       0;
            A.m(3)    0     A.m(1)    0       0       0       0       0;
            A.m(4)    0       0     A.m(1)    0       0       0       0;
            A.m(5) -A.m(3)  A.m(2)    0     A.m(1)    0       0       0;
            A.m(6)    0    -A.m(4)  A.m(3)    0     A.m(1)    0       0;
            A.m(7)  A.m(4)    0    -A.m(2)    0       0     A.m(1)    0;
            A.m(8)  A.m(6)  A.m(7)  A.m(5)  A.m(4)  A.m(2)  A.m(3)  A.m(1)];
            E = expm(M);
            r = PGA(E(1:8,1));
        end

        function r = bilinear_(A, B)
            %bilinear(A, B): Computes the bilinear form of two GA objects.

            S = GAsignature;
            N = [B.m(1);S(1)*B.m(2);S(2)*B.m(3);S(3)*B.m(4);S(1)*S(2)*B.m(5);S(2)*S(3)*B.m(6);S(3)*S(1)*B.m(7);S(1)*S(2)*S(3)*B.m(8)];
            r = PGA((A.m')*N);
        end

        function r = connection(e,A,B)
            %connection(e,A,B): Compute the translational moment from A to B.
            %  Both A and B must be in homogeneous form.
            %  This is what you need to add to A to make its meet with B non-trivial.
            e = PGA(e);
            A = PGA(A);
            B = PGA(B);

            tA = inner(e,A);
            tB = inner(e,B);
            ejAB = e^(join(tA,tB)*meet(tA,tB));
            Ato0 = PGAZ((A-contraction(A,ejAB)/ejAB)/inner(e,A));
            Bto0 = PGAZ((B-contraction(B,ejAB)/ejAB)/inner(e,B));

            r = PGA.returnPGA(PGAZ(PGA(Bto0-Ato0))^tA);
        end

        % ***** Draw functions *****

        function r = PGAtext(v,s,c)
            %r = PGAtext(v,s,c): draw string s at tip of vector v in color c (optional).
            %  r is the text handle.

            % TODO: test for other nargin values
            if nargin==2
                r = text(v.m(2),v.m(3),v.m(4),s);
            else
                r = text(v.m(2),v.m(3),v.m(4),s,'Color',c);
            end
        end

        function r = geoall(A,B)
            %geoall(A,B): Compute the geometric relationships between blades A and B
            % Return a structure with the following fields
            %   obj1 = A
            %   obj2 = B
            %   comp = orthogonal complement of A in B
            %   proj = projection of A onto B
            %   rej  = rejection of A by B
            %   meet = meet (i.e. intersection) of A and B
            %   join = join (i.e. union) of A and B
            
            if grade(PGA(A)) == -1 & grade(PGA(B)) == -1
                error('A and B should be pure blades')
            end
            r.obj1 = A;
            r.obj2 = B;
            r.comp = contraction(A,B);  
            r.proj = r.comp/B;
            r.rej = A - r.proj;
            r.meet = meet(A,B);
            r.join = join(A,B);
        end

        % ***** Public Drawing Methods *****
        % ZNOTE: not sure about whether it makes sense to return bivectors, trivectors, as part of the drawing method

        function b = DrawBivector(A,B,c)
            % DrawBivector(A,B,c) : draw the bivector A^B.  
            %  A and B must be vectors.  c is an optional color argument.
            
            if nargin == 2
                c = 'y';
            end

            % ZNOTE: TODO: Think about whether we should wrap PGAZ with PGA.
            % Also, should PGAZ be private?
            A = PGA(PGAZ(A));
            B = PGA(PGAZ(B));

            if ~PGAisa(A,'vector') || ~PGAisa(B,'vector')
                % ZNOTE: TODO: Should we be calling display here instead?
                A
                B
                error('DrawBivector: A and B must both be vectors.');
            end
            % DrawBivector: draw the parallelogram bounded by A,B.
            X = [0 A.m(2) A.m(2)+B.m(2) B.m(2)];
            Y = [0 A.m(3) A.m(3)+B.m(3) B.m(3)];
            Z = [0 A.m(4) A.m(4)+B.m(4) B.m(4)];
            PGA.biarrow(B, A, 'g');
            patch(X,Y,Z,c);
            draw(A);
            b = A^B;
        end

        function t = DrawTrivector(A,B,C)
            % DrawTrivector(A,B,C): draw the parallelepiped bounded by A,B,C.  
            %  First call DrawBivector(A,B) and then add the rest of the arrows.
            draw(A,'b');
            draw(B,'g');
            draw(C,'m');
            X = [C.m(2) C.m(2)+A.m(2) C.m(2)+A.m(2)+B.m(2) C.m(2)+B.m(2) C.m(2)];
            Y = [C.m(3) C.m(3)+A.m(3) C.m(3)+A.m(3)+B.m(3) C.m(3)+B.m(3) C.m(3)];
            Z = [C.m(4) C.m(4)+A.m(4) C.m(4)+A.m(4)+B.m(4) C.m(4)+B.m(4) C.m(4)];
            plot3(X,Y,Z,'r');
            X = [C.m(2)+A.m(2) A.m(2) A.m(2)+B.m(2) C.m(2)+A.m(2)+B.m(2)];
            Y = [C.m(3)+A.m(3) A.m(3) A.m(3)+B.m(3) C.m(3)+A.m(3)+B.m(3)];
            Z = [C.m(4)+A.m(4) A.m(4) A.m(4)+B.m(4) C.m(4)+A.m(4)+B.m(4)];
            plot3(X,Y,Z,'r');
            X = [A.m(2)+B.m(2)  B.m(2) C.m(2)+B.m(2)]';
            Y = [A.m(3)+B.m(3) B.m(3) C.m(3)+B.m(3)]';
            Z = [A.m(4)+B.m(4) B.m(4) C.m(4)+B.m(4)]';
            plot3(X,Y,Z,'r');
            axis equal;
            t = A^B^C;
        end

        function draw(s, O, c)
            %draw(s, O, c) - draw a multivector. 
            %  s: the multivector to be drawn.
            %  O: an offset by which to displace the multivector before drawing (optional)
            %  c: a color in which to draw the multivector (optional).
            %
            %  For bivectors, can set shape using 'GAbvShape'

            s = PGA(s);
            

            set(gcf, 'Render', GArender)
            
            A = PGAZ(s);
            
            %ZNOTE: I don't like this way this is done. Find a better way.

            % if you assign to nargin, then Matlab doesn't set it
            % We will set na to be 2 (no color) or 3 (color)
            na = nargin;	
            
            if nargin == 1
                O = PGA([0;0;0;0;0;0;0;0]);
                na = 2;
            elseif nargin == 2
                if isa(O,'char')
                   c = O;
                   O = PGA([0;0;0;0;0;0;0;0]);
                   na = 3;
                else
                    O = PGA(O);
                end
            else
                O = PGA(O);
            end
            
            if isa(A,'double') 
                if nargin == 2
                 title(['scalar = ' num2str(A)]);
                else
                 title(['scalar = ' num2str(A)], 'Color', c);
                end
                hold on;
            elseif PGAisa(A,'double')
                if nargin == 2
                 title(['scalar = ' num2str(A.m(1))]);
                else
                 title(['scalar = ' num2str(A.m(1))], 'Color', c);
                end
                hold on;
            elseif PGAisa(A,'vector')
                 if na == 2
                    PGA.arrow(A,O);
                 else
                    PGA.arrow(A,O,c);
                 end
                 hold on;
            elseif PGAisa(A,'bivector')
                 if na == 2
                    PGA.drawBivector(A,O);
                 else
                    PGA.drawBivector(A,O,c);
                 end
            elseif PGAisa(A,'trivector')
                 if na == 2
                    PGA.drawTrivector(A,O);
                 else
                    PGA.drawTrivector(A,O,c);
                 end
            else
                 % Draw a multivector as the sum of its parts
                 if sum(abs(A.m(2:4))) > 0
                      M = PGA([0; A.m(2:4); 0; 0; 0; 0]);
                  if na == 3
                    draw(M,O,c);
                      else
                    draw(M,O);
                  end
                 end
                 if sum(abs(A.m(5:7))) > 0
                      M = PGA([0; 0; 0; 0; A.m(5:7); 0]);
                  if na == 3
                    draw(M,O,c);
                      else
                    draw(M,O);
                  end
                 end
                 if A.m(8) ~= 0
                      M = PGA([0; 0; 0; 0; 0; 0; 0; A.m(8)]);
                  if na == 3
                    draw(M,O,c);
                      else
                    draw(M,O);
                  end
                 end
                 if A.m(1) ~= 0
                  if na == 3
                    draw(A.m(1),O,c);
                      else
                    draw(A.m(1),O);
                  end
                 end
                 return;	% Return to avoid 2 'axis' calls
            end
            % We want equal axes to avoid distortion.  However, we need
            % the axis('auto') before hand to ensure that everything
            % appears on the screen.
            axis('auto');
            axis('equal');
            axis('tight');
        end     
        
        function r = DrawOuter(A, B)
            %DrawOuter(A,B): Draw A^B

            r = PGA.DrawProduct(A, B, 'wedge');
        end    
        
        function r = DrawInner(A, B)
            %DrawInner(A,B): draw the inner product of A,B.

            %ZNOTE: Shouldn't DrawProduct be lower case? What is the convention here?
            r = PGA.DrawProduct(A, B, 'inner');
        end

        function r = DrawGP(A, B)
            %DrawGP(A,B): Draw the geometric product A*B.
            
            r = PGA.DrawProduct(A, B, 'geometric');
        end

        % ZNOTE: unclear why this should be a method within PGA itself.
        function plot3(A, c)
            % plot3(A, c): plot a 3d line connecting end points of A in color c
            
            if nargin == 1
                c = 'b';
            end
            for i = 1:length(A)
                p = A{i};
                x(i) = p.m(1);
                y(i) = p.m(2);
                z(i) = p.m(3);
            end
            plot3(x,y,z,c);
        end
    end






















    methods (Access = public)
        % Constructor.
        % if no argument is supplied, creates the zero multivector
        % if a multivector is supplied, returns that multivector
        % if two integers are supplied, returns 
        % TODO: finish the above comment
        function obj = PGA(m0, m1, m2, m3, m4)
            if nargin == 0
                obj = PGA(0);
            elseif nargin == 5
                obj.m = [m0; m1(1); m1(2); m1(3); m1(4); m2(1); m2(2); m2(3); m2(4); m2(5); m2(6); m2(7); m2(8); m2(9); m2(10); m2(11); m2(12); m3(1); m3(2); m3(3); m3(4); m4];
            elseif isa(m0, 'PGA')
                obj = m0;
            elseif isa(m0, 'double')
                if size(m0, 1) == 1 & size(m0, 2) == 1
                    % User has provided a scalar
                    % TODO: constructor here is outdated
                    obj.m = [m0; 0; 0; 0; 0; 0; 0; 0];
                elseif size(m0, 1) == 1 & size(m0, 2) == 8
                    % User has provided a column vector
                    obj.m = m0';
                elseif size(m0, 1) == 8 & size(m0, 2) == 1
                    % User has prodived a row vector
                    obj.m = m0;
                else
                    error('Bad GA argument array size.\nExpected size is either 1x1, 8x1 or 1x8.\nCurrent size is: %dx%d', size(m, 1), size(m, 2))
                end
            else
                error('Bad PGA argument: must be an array of doubles. Class is currently: %s.', class(m))
            end
        end

        % Returns the matrix for the PGA object. For debugging purposes.
        % TODO: add to change notes: method m in GA is now called matrix in PGA.
        function r = matrix(self)
            r = self.m;
        end

        function r = PGAisa(A, t)
            % TODO: check if t is the correct type.
            r = PGA(A).PGAisa_(t);
        end

        function r = double(A)
            r = double_(PGA(A));
        end

        % Overloaded methods. ZNOTE: rename this section

        function r = plus(A, B)
            r = PGA.returnPGA(plus_(PGA(A), PGA(B)));
        end

        function r = minus(A, B)
            r = PGA.returnPGA(minus_(PGA(A), PGA(B)));
        end

        % ZNOTE: does not have a private counterpart.
        function r = uminus(A)
            r = PGA.returnPGA(-PGA(A).m);
        end

        function r = mpower(A, B)
            r = outer(A, B);
        end

        function r = outer(A, B)
            r = PGA.returnPGA(outer_(PGA(A), PGA(B)));
        end

        function r = inner(A, B)
            r = PGA.returnPGA(inner_(PGA(A), PGA(B)));
        end

        function r = contraction(A, B)
            r = PGA.returnPGA(contraction_(PGA(A), PGA(B)));
        end

        function r = innerS(A, B)
            r = PGA.returnPGA(innerS_(PGA(A), PGA(B)));
        end

        function r = innerH(A, B)
            r = PGA.returnPGA(innerH_(PGA(A), PGA(B)));
        end

        % ***** Printing Methods *****

        function display(p)
            disp(' ');
            disp([inputname(1),' = '])
            disp(' ');
            disp([char(p)])
            disp(' ');
        end
        
        % ***** Equality and Inequality *****
        function r = eq(A, B)
            r = eq_(PGA(A), PGA(B));
        end

        function r = eeq(A, B)
            r = eeq_(PGA(A), PGA(B));
        end

        function r = ne(A, B)
            r = ~eq_(PGA(A), PGA(B));
        end

        function r = product(A, B)
            r = PGA.returnPGA(product_(PGA(A), PGA(B)));
        end

        function r = mtimes(A, B)
            r = product(A, B);
        end

        function r = inverse(A)
            r = PGA.returnPGA(inverse_(PGA(A)));
        end

        function r = mrdivide(A, B)
            r = PGA.returnPGA(PGA(A) * inverse_(PGA(B)));
        end

        function r = bilinear(A, B)
            %bilinear(A, B): Computes the bilinear form of two GA objects.
            r = PGA.returnPGA(bilinear_(PGA(A), PGA(B)));
        end

        function r = wexp(A)
            r = PGA.returnPGA(wexp_(PGA(A)));
        end

        function r = gexp(A)
            r = PGA.returnPGA(gexp_(PGA(A)));
        end

        function r = norm(A)
            r = norm_(PGA(A));
        end

        function r = dual(A)
            r = PGA.returnPGA(dual_(PGA(A)));
        end

        function r = reverse(A)
            r = PGA.returnPGA(reverse_(PGA(A)));
        end

        function r = meet(A, B)
            r = PGA.returnPGA(meet_(PGA(A), PGA(B)));
        end

        function r = join(A, B)
            r = PGA.returnPGA(join_(PGA(A), PGA(B)));
        end

        function r = sLog(spinor)
            % TODO: Spinor verification
            r = sLog_(spinor);
        end

        function r = conjugate(A)
            r = PGA.returnPGA(conjugate_(PGA(A)));
        end

        function r = gradeinvolution(A)
            r = PGA.returnPGA(gradeinvolution(PGA(A)));
        end

        function r = grade(A, n)
            % TODO: ensure correct type of n
            r = PGA.returnPGA(grade_(PGA(A), n));
        end

        function r = isGrade(A)
            r = isGrade_(PGA(A));
        end
    end
end
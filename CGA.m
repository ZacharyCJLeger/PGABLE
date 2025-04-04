classdef CGA < GA
    %CGA - A child class of GA for elements of Conformal Geometric Algebra.
    %   Elements
    %      Basic elements include 1, e, e1, e2, e3, eb, no, ni, e12, e23, e31, 
    %      e123
    %      Additionally, we have e21 = -e12, e32 = -e23, e13 = -e31
    %      We also have method for creating CGA points,gapoint(x, y, z), which creates a
    %      CGA point with coordinates (x, y, z). We also have origin() = gapoint(0, 0, 0).
    %
    %   Operations
    %      You can use these special characters for these basic operations:
    %         • +  for addition               also: plus(A, B)
    %         • -  for subtraction            also: minus(A, B)
    %         • *  for the geometric product  also: product(A, B)
    %         • /  for division               also: divide(A, B)
    %         • ^  for the outer product      also: outer(A, B)
    %         • .* for the inner product      also: inner(A, B)
    %         • == for equality               also: eq(A, B)
    %         • ~= for inequality             also: neq(A, B)
    %      Additonally, there are basic operations:
    %         • meet(A, B)                    compute the meet of two multivectors
    %         • join(A, B)                    compute the join of two multivectors
    %         • inverse(A)                    compute the inverse
    %           (Note: the inverse may not always exist in CGA)
    %         • gradeinvolution(A)            compute the grade involution
    %         • conjugate(A)                  compute the conjugate
    %         • reverse(A)                    compute the reverse
    %         • norm(A)                       compute the norm
    %         • vnorm(A)                      compute the vanishing norm
    %         • normalize(A)                  normalize the multivector
    %         • dual(A)                       compute the dual
    %         • getx(A)                       get the x coordinate of a CGA point
    %         • gety(A)                       get the y coordinate of a CGA point
    %         • getz(A)                       get the z coordinate of a CGA point
    %         • zeroepsilons(A)               zero-out epsilons (small errors)
    %         • draw(A)                       draw the multivector
    %         (See also GAScene for more information on draw calls)
    %         • grade(A, g)                   select the grade-g component of a multivector
    %         • isgrade(A, g)                 determine if a multivector is of grade g
    %      There are functions for constructing some objects directly:
    %         •gapoint(x,y,z)                 construct a CGA point
    %         •galine(l,p)                    construct a line with direction
    %         •galine(l1,l2,l3, p1,p2,p3)     with direction l through point p
    %      There are also more advanced operations:
    %         • sqrt(A)                       compute the square root
    %         • glog(A)                       compute the geometric log
    %         • gexp(A)                       compute the geometric exponent
    %
    %   See also GA, OGA, PGA, GAScene.

    % PGABLE, Copyright (c) 2025, University of Waterloo
    % Copying, use and development for non-commercial purposes permitted.
    %          All rights for commercial use reserved; for more information
    %          contact Stephen Mann (smann@uwaterloo.ca)
    %
    %          This software is unsupported.
    
    properties (Access = private)
        % A 1x32 matrix of real numbers corresponding to the coefficients of entries 1, e0, e1, ..., e01, ..., e0123. 
        m
    end

    %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%
    %           Settings           %
    %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%

    methods (Static = true)

        function settings()
            %SETTINGS - Displays the current configuration settings for CGA in PGABLE.
            %   To retrieve a particular setting, run CGA.[setting].
            %   For example, to retrieve the value of increasing_order, run
            %   PGA.increasing_order.
            %   To change the value of a particular setting, run CGA.[setting]([value]).
            %   For example, to set the value of increasing_order to true, run
            %   CGA.increasing_order(true).
            %   For more information on a particular setting, run help CGA.[setting].
            %
            %   See also GA.settings.

            [S0, S1, S2, S3, S4] = CGA.signature();

            disp("   ~~~~~~~~~~ CGA Settings ~~~~~~~~~~")
            disp("   signature:        no*no = " + S0 + ", e1*e1 = " + S1 + ", e2*e2 = " + S2 + ", e3*e3 = " + S3 + ", ni*ni = " + S4)
            disp("   increasing_order: " + CGA.increasing_order())
            disp("   ~~~~~~~~~~ CGA Point Settings ~~~~~~~~~~")
            disp("   pointsize:        " + CGA.pointsize())
        end
        
        function [S0, S1, S2, S3, S4] = signature(sign0, sign1, sign2, sign3, sign4)
            %SIGNATURE - Set/retrieve the current signature of the model.
            %   This setting is NOT recommended for beginners.
            %   If no arguments are provided, the signatures for no, e1, e2, e3, ni are returned
            %   as a vector [S0, S1, S2, S3, S4].
            %   If 5 arguments are provided, the inputs sign0, sign1, sign2, sign3, sign4 correspond
            %   to the signatures of no, e1, e2, e3, ni respectively.

            persistent signature0;
            persistent signature1;
            persistent signature2;
            persistent signature3;
            persistent signature4;

            if isempty(signature0)
                signature0 = 0;
                signature1 = 1;
                signature2 = 1;
                signature3 = 1;
                signature4 = 0;
            end

            if nargin == 5
               signature0 = sign0;
               signature1 = sign1;
               signature2 = sign2;
               signature3 = sign3;
               signature4 = sign4;
            end

            S0 = signature0;
            S1 = signature1;
            S2 = signature2;
            S3 = signature3;
            S4 = signature4;
        end

        function val = pointsize(newval, surpress_output)
            %POINTSIZE - Set/retreive the POINTSIZE setting.
            %   The POINTSIZE setting is a double indicated the radius of the octahedron
            %   representing a point.
            %   If no argument is provided, POINTSIZE returns the current value of the
            %   POINTSIZE setting.
            
            arguments
                newval = [];
                surpress_output = false;
            end

            persistent currentval;
            
            if isempty(currentval)
                currentval = 0.1;
            end

            if isempty(newval)
                % User is trying to retrieve the current value
                val = currentval;
            else
                % User is trying to set the value
                if isnumeric(newval)
                    currentval = newval;
                    if ~surpress_output
                        disp("   pointsize set to " + currentval)
                    end
                else
                    error('pointsize must have a numeric value.')
                end
            end 
        end

        function val = increasing_order(newval, surpress_output)
            %INCREASING_ORDER - Set/retrieve the INCREASING_ORDER setting.
            %   The INCREASING setting is either true or false.
            %   When set to true, PGA elements are represented by the basis:
            %   1, e0, e1, e2, e3, e01, e02, e03, e12, e13, e23, e012, e013, e023, e123, e0123
            %   When set to false, PGA elements are represented by the basis:
            %   1, e0, e1, e2, e3, e01, e02, e03, e23, e31, e12, e032, e013, e021, e123, e0123
            %   If no argument is provided, INCREASING_ORDER returns the current value of the
            %   INCREASING_ORDER setting.

            arguments
                newval = [];
                surpress_output = false;
            end

            persistent currentval;
            
            % By default the increasing_order setting is set to false
            if isempty(currentval)
                currentval = false;
            end

            if isempty(newval)
                % User is trying to retrieve the current value
                val = currentval;
            else
                % User is trying to set the value
                if islogical(newval)
                    currentval = newval;
                    if ~surpress_output
                        disp("   increasing_order set to " + currentval)
                    end
                else
                    error('increasing_order must have value true or false.')
                end
            end
        end

        %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%
        %  Dynamic Drawing Functions   %
        %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%

        function h = drawvanishingpoint(vp, varargin)
            hold on


            vpx = vp.getx();
            vpy = vp.gety();
            vpz = vp.getz();
            
            ax = gca;

            xrange = ax.XLim;
            yrange = ax.YLim;
            zrange = ax.ZLim;

            xwidth = xrange(2) - xrange(1);
            ywidth = yrange(2) - yrange(1);
            zwidth = zrange(2) - zrange(1);

            propx = abs(vpx/xwidth);
            propy = abs(vpy/ywidth);
            propz = abs(vpz/zwidth);

            xaverage = (xrange(1) + xrange(2))/2;
            yaverage = (yrange(1) + yrange(2))/2;
            zaverage = (zrange(1) + zrange(2))/2;

            center = gapoint(xaverage, yaverage, zaverage, CGA);

            
            [~, argmax] = max([propx, propy, propz]);

            switch argmax
                case 1
                    if vp.getx() > 0
                        plane = -xrange(2)*e0(PGA) + e1(PGA);
                    else
                        plane = -xrange(1)*e0(PGA) + e1(PGA);
                    end
                case 2
                    if vp.gety() > 0
                        plane = -yrange(2)*e0(PGA) + e2(PGA);
                    else
                        plane = -yrange(1)*e0(PGA) + e2(PGA);
                    end
                case 3
                    if vp.getz() > 0
                        plane = -zrange(2)*e0(PGA) + e3(PGA);
                    else
                        plane = -zrange(1)*e0(PGA) + e3(PGA);
                    end
            end

            line = join(vp, center);
            line = normalize(line);
            p = line^plane;

            % Phi is percent away the arrows are from the edge of the bounding box
            phi = 0.3;
            aphi = 1 - phi;
            h = plot3([p.getx(), aphi*p.getx() + phi*xaverage], ...
                      [p.gety(), aphi*p.gety() + phi*yaverage], ...
                      [p.getz(), aphi*p.getz() + phi*zaverage], 'k');


            h = [h, PGABLEDraw.drawstaronplane(vp, 1, varargin{:})];
            h = [h, PGABLEDraw.drawstaronplane(vp, 2, varargin{:})];
            h = [h, PGABLEDraw.drawstaronplane(vp, 3, varargin{:})];
        end

        function h = drawvanishingline(vl, varargin)
            %DRAWVANISHINGLINE - Draws a single instance of a vanishing line.
            %   This function is NOT intended for a user to draw a vanishing line to the
            %   scene. To draw a vanishing line, run "draw(vanishing_line)".

            h = [];

            ax = gca;

            xrange = ax.XLim;
            yrange = ax.YLim;
            zrange = ax.ZLim;

            xaverage = (xrange(1) + xrange(2))/2;
            yaverage = (yrange(1) + yrange(2))/2;
            zaverage = (zrange(1) + zrange(2))/2;
            
            % TODO: move this somewhere else
            hold on
            
            % TODO: get x, y, z of line's normal.
            vlplane = normalize(hd(vl)/I3(PGA));

            n = normalize(hd(vl)/I3(PGA));
            
            p0 = xaverage*e1(PGA) + yaverage*e2(PGA) + zaverage*e3(PGA);

            vlplane = n - (n.*p0)*e0;
            


            % The 6 planes of the bounding box
            xp = -xrange(2)*e0(PGA) + e1(PGA);
            xn = -xrange(1)*e0(PGA) + e1(PGA);
            yp = -yrange(2)*e0(PGA) + e2(PGA);
            yn = -yrange(1)*e0(PGA) + e2(PGA);
            zp = -zrange(2)*e0(PGA) + e3(PGA);
            zn = -zrange(1)*e0(PGA) + e3(PGA);
            

            points = [];

            isin = @(p)PGABLEDraw.isinboundingbox(xrange, yrange, zrange, p);

            if isin(xp^yp^vlplane)
                % We are in the scenario where we know the line goes through:
                %   - The x planes and the y planes
                % However, it may still go through the z planes.
                % If it does, it must hit the following point:

                if isin(yp^zp^vlplane)
                    % Now we know we go through all six faces in this particular order.
                    points = {xp^yp^vlplane, yp^zp^vlplane, zp^xn^vlplane, xn^yn^vlplane, yn^zn^vlplane, zn^xp^vlplane, xp^yp^vlplane};
                elseif isin(yp^zn^vlplane)
                    points = {xp^yp^vlplane, yp^zn^vlplane, zn^xn^vlplane, xn^yn^vlplane, yn^zp^vlplane, zp^xp^vlplane, xp^yp^vlplane};
                else
                    % It doesn't go through all 6 faces. So we get this
                    points = {xp^yp^vlplane, yp^xn^vlplane, xn^yn^vlplane, yn^xp^vlplane, xp^yp^vlplane};
                end
            elseif isin(xp^zp^vlplane)
                if isin(xp^yn^vlplane)
                    points = {zp^xp^vlplane, xp^yn^vlplane, yn^zn^vlplane, zn^xn^vlplane, xn^yp^vlplane, yp^zp^vlplane, zp^xp^vlplane};
                else
                    points = {xp^zp^vlplane, zp^xn^vlplane, xn^zn^vlplane, zn^xp^vlplane, xp^zp^vlplane};
                end
            else
                if isin(xp^zn^vlplane)
                    points = {zn^xp^vlplane, xp^yn^vlplane, yn^zp^vlplane, zp^xn^vlplane, xn^yp^vlplane, yp^zn^vlplane, zn^xp^vlplane};
                else
                    points = {yp^zp^vlplane, zp^yn^vlplane, yn^zn^vlplane, zn^yp^vlplane, yp^zp^vlplane};
                end
            end

            points = arrayfun(@(p)PGABLEDraw.boundingboxclip(xrange, yrange, zrange, p), points);

            h = PGABLEDraw.plotline(points, varargin{:});

            
            for i=1:(length(points)-1)
                point_1 = points{i};
                point_2 = points{i+1};
                % TODO: This is very bad. Very bad. Getting imaginary parts. Draws the points correctly though.
                ap_move = gexp(glog(point_1/point_2)*0.5);
                ap = sqrt(ap_move)*point_2/sqrt(ap_move);
                % Now we want to draw the line from AP to AP moved in the direction where things will go
                % Then, we point with an arrow type thing.
                ap_short_move = gexp(glog(point_1/point_2)*0.45);
                ap_long_move = gexp(glog(point_1/point_2)*0.55);
                ap_short = sqrt(ap_short_move)*point_2/sqrt(ap_short_move);
                ap_long = sqrt(ap_long_move)*point_2/sqrt(ap_long_move);

                R = gexp(-0.05*vl/2);
                ap_tip = R*ap/R;

                arrow_points = {ap_short, ap_tip, ap_long};
                arrow_points = arrayfun(@(p)PGABLEDraw.boundingboxclip(xrange, yrange, zrange, p), arrow_points);
                
                c = PGABLEDraw.extractcolor(varargin{:});
                h = [h PGABLEDraw.patch(arrow_points, 'EdgeColor', c, 'FaceColor', c)];
            end
            
        end
    end

    %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%
    %      Protected methods       %
    %         (non-static)         %
    %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%
    
    methods (Access = protected)

        function [scalar_nz, vector_nz, bivector_nz, trivector_nz, fourvector_nz, fivevector_nz] = gradestatus_(A)
            nzm = A.m ~= 0;
            % The following variables hold true if there exists a non-zero entry
            % in that category
            scalar_nz = nzm(1) ~= 0;
            vector_nz = sum(nzm(2:6)) ~= 0;
            bivector_nz = sum(nzm(7:16)) ~= 0;
            trivector_nz = sum(nzm(17:26)) ~= 0;
            fourvector_nz = sum(nzm(27:31)) ~= 0;
            fivevector_nz = nzm(32) ~= 0;
        end

        function b = GAisa_(A, t)
            [scalar_nz, vector_nz, bivector_nz, trivector_nz, fourvector_nz, fivevector_nz] = gradestatus_(A);

            [scal, ...
            EO, E1, E2, E3, EI, ...
            EO1, EO2, EO3, EOI, E12, E13, E1I, E23, E2I, E3I, ...
            EO12, EO13, EO1I, EO23, EO2I, EO3I, E123, E12I, E13I, E23I, ...
            EO123, EO12I, EO13I, EO23I, E123I, ...
            EO123I] = decompose_(A);
	     
            if strcmp(t, 'point') || strcmp(t,'sphere') || strcmp(t, 'plane')
                if isgrade(A,1)
                    if abs(A.m(2)) < GA.epsilon_tolerance()
                        if strcmp(t, 'plane') 
                            b = true;
                        else
                            b = false;
                        end
                        return;
                    end

	                ptc = A.m(2:6)/A.m(2);
		            v = ptc(2:4);
		            if strcmp(t, 'point')
                        %norm(v)
                        %norm(v)*norm(v)/2 - ptc(5)
                        %norm(v)*GA.epsilon_tolerance()
                        % We need the scaling of epsilon just because that's how the numerics work
                        % The '=' is to handle the origin
                        if abs(norm(v)*norm(v)/2 - ptc(5)) <= 10*norm(v)*GA.epsilon_tolerance() 
                            b = true;
                        else
                            b = false;
                        end
                    else
                        b = true;
                    end
		            return;
	            elseif isgrade(A,4)
		            A = dual(A);
                    if strcmp(t,'sphere')
                        if abs(A.m(2)) > GA.epsilon_tolerance()
                            b = true;
                        else
                            b = false;
                        end
                        return;
		            elseif strcmp(t,'plane')
                        if abs(A.m(2)) < GA.epsilon_tolerance()
                                b = true;
                        else
                            b = false;
                        end
		                return;
                    end
                else
	                b = false;
		        return;
	            end
	        end

            if strcmp(t, 'circle') || strcmp(t,'line')
                if isgrade(A,2)
                    if  abs(E3I)< GA.epsilon_tolerance()  &&  abs(E2I)< GA.epsilon_tolerance()  &&  abs(E1I)< GA.epsilon_tolerance()
                        if strcmp(t,'line')
                            b = true;
                        else
                            b = false;
                        end
                        return
                    else
                        if strcmp(t,'circle')
                            b = true;
                        else
                            b = false;
                        end
                        return
                    end
                elseif isgrade(A,3)
                    if  abs(EO12) < GA.epsilon_tolerance() && abs(EO13) < GA.epsilon_tolerance() && abs(EO23) < GA.epsilon_tolerance()
                        if strcmp(t,'line')
                            b = true;
                        else
                            b = false;
                        end
                        return
		            else
                        if strcmp(t,'circle')
                            b = true;
                        else
                            b = false;
                        end
		            return
		            end
                else
                    b = false;
                return;
                end
            end
	    
            if strcmp(t, 'double') || strcmp(t, 'scalar') 
                b = ~(             vector_nz || bivector_nz || trivector_nz || fourvector_nz || fivevector_nz);
            elseif strcmp(t,'vector') 
                b = ~(scalar_nz ||              bivector_nz || trivector_nz || fourvector_nz || fivevector_nz);
            elseif strcmp(t,'bivector') || strcmp(t,'line')
                b = ~(scalar_nz || vector_nz ||                trivector_nz || fourvector_nz || fivevector_nz);
            elseif strcmp(t,'trivector') 
                b = ~(scalar_nz || vector_nz || bivector_nz ||                 fourvector_nz || fivevector_nz);
            elseif strcmp(t,'4vector') || strcmp(t,'fourvector') || strcmp(t,'quadvector') 
                b = ~(scalar_nz || vector_nz || bivector_nz || trivector_nz                  || fivevector_nz );
            elseif strcmp(t,'5vector') || strcmp(t,'fivevector') || strcmp(t,'quintvector') || strcmp(t,'pseudoscalar')
                b = ~(scalar_nz || vector_nz || bivector_nz || trivector_nz || fourvector_nz                  );
            elseif strcmp(t,'multivector')
                b = sum([scalar_nz vector_nz bivector_nz trivector_nz fourvector_nz fivevector]) > 1;
            elseif strcmp(t,'plane')
  	        % this is wrong
                b = ~(scalar_nz ||              bivector_nz || trivector_nz || fourvector_nz || fivevector_nz);
            else
                b = false;
            end 
        end

        function r = double_(A)
            if GAisa_(A, 'scalar')
                r = A.m(1);
            else
                error('Can only convert a scalar CGA object to a double. Object is %s.', char(A));
            end
        end

        % ***** Functions for adding and subtracting CGA objects *****

        function R = plus_(A, B)
            R = CGA(A.m + B.m);
        end

        function R = minus_(A, B)
            R = CGA(A.m - B.m);
        end

        function R = uminus_(A)
            R = CGA(-A.m);
        end


        function [SO,S1,S2,S3,SI, ...
                SO1, SO2, SO3, SOI, S12, S13, S1I, S23, S2I, S3I, ...
                SO12, SO13, SO1I, SO23, SO2I, S123, S12I, S23I,...
                SO123, SO12I, SO13I, SO23I, S123I, SO123I] = computeSCoef()
            [SO, S1, S2, S3, SI] = CGA.signature();
            SO1 = SO*S1;
            SO2 = SO*S2;
            SO3 = SO*S3;
            SOI = SO*SI;
            S12 = S1*S2;
            S13 = S1*S3;
            S1I = S1*SI;
            S23 = S2*S3;
            S2I = S2*SI;
            S3I = S3*SI;
            
            SO12 = SO1*S2;
            SO13 = SO1*S3;
            SO1I = SO1*SI;
            SO23 = SO2*S3;
            SO2I = SO2*SI;
            S123 = S12*S3;
            S12I = S12*SI;
            S23I = S23*SI;

            SO123 = SO12*S3; 
            SO12I = SO12*SI; 
            SO13I = SO13*SI; 
            SO23I = SO23*SI; 
            S123I = S123*SI;
            
            SO123I = SO123*SI; 
        end
        % ***** Geometric Product Stuff *****

        function R = product_(A, B)
            [scal, ...
            EO, E1, E2, E3, EI, ...
            EO1, EO2, EO3, EOI, E12, E13, E1I, E23, E2I, E3I, ...
            EO12, EO13, EO1I, EO23, EO2I, EO3I, E123, E12I, E13I, E23I, ...
            EO123, EO12I, EO13I, EO23I, E123I, ...
            EO123I] = decompose_(A);

% 	    [SO,S1,S2,S3,SI, SO1, SO2, SO3, SOI, S12, S13, S1I, S23, S2I, S3I, SO12, SO13, SO1I, SO23, SO2I, S123, S12I, S23I, SO123, SO12I, SO13I, SO23I, S123I, SO123I] = computeSCoef();

        
	Am = [
        scal          -EI           E1           E2           E3          -EO         -E1I         -E2I         -E3I          EOI         -E12         -E13         -EO1         -E23         -EO2         -EO3         E12I         E13I         EO1I         E23I         EO2I         EO3I        -E123         EO12         EO13         EO23        E123I       -EO12I       -EO13I       -EO23I        EO123      -EO123I ;
          EO  scal + -EOI          EO1          EO2          EO3            0  -E1 + -EO1I  -E2 + -EO2I  -E3 + -EO3I           EO        -EO12        -EO13            0        -EO23            0            0  -E12 + EO12I  -E13 + EO13I         -EO1  -E23 + EO23I         -EO2         -EO3       -EO123            0            0            0  E123 + EO123I        -EO12        -EO13        -EO23            0        EO123 ;
          E1         -E1I         scal          E12          E13          EO1          -EI        -E12I        -E13I        -EO1I          -E2          -E3           EO        -E123         EO12         EO13          E2I          E3I         -EOI        E123I       -EO12I       -EO13I         -E23         -EO2         -EO3       -EO123         E23I         EO2I         EO3I       EO123I        -EO23        EO23I ;
          E2         -E2I         -E12         scal          E23          EO2         E12I          -EI        -E23I        -EO2I           E1         E123        -EO12          -E3           EO         EO23         -E1I       -E123I        EO12I          E3I         -EOI       -EO23I          E13          EO1        EO123         -EO3        -E13I        -EO1I      -EO123I         EO3I         EO13       -EO13I ;
          E3         -E3I         -E13         -E23         scal          EO3         E13I         E23I          -EI        -EO3I        -E123           E1        -EO13           E2        -EO23           EO        E123I         -E1I        EO13I         -E2I        EO23I         -EOI         -E12       -EO123          EO1          EO2         E12I       EO123I        -EO1I        -EO2I        -EO12        EO12I ;
          EI            0         -E1I         -E2I         -E3I   scal + EOI            0            0            0          -EI        -E12I        -E13I   E1 + -EO1I        -E23I   E2 + -EO2I   E3 + -EO3I            0            0         -E1I            0         -E2I         -E3I        E123I  -E12 + -EO12I  -E13 + -EO13I  -E23 + -EO23I            0         E12I         E13I         E23I  -E123 + EO123I        E123I ;
         EO1  -E1 + -EO1I           EO         EO12         EO13            0  scal + -EOI  E12 + -EO12I  E13 + -EO13I          EO1         -EO2         -EO3            0       -EO123            0            0    E2 + EO2I    E3 + EO3I          -EO  E123 + EO123I        -EO12        -EO13        -EO23            0            0            0  -E23 + EO23I         -EO2         -EO3       -EO123            0         EO23 ;
         EO2  -E2 + -EO2I        -EO12           EO         EO23            0  -E12 + EO12I  scal + -EOI  E23 + -EO23I          EO2          EO1        EO123            0         -EO3            0            0  -E1 + -EO1I  -E123 + -EO123I         EO12    E3 + EO3I          -EO        -EO23         EO13            0            0            0  E13 + -EO13I          EO1        EO123         -EO3            0        -EO13 ;
         EO3  -E3 + -EO3I        -EO13        -EO23           EO            0  -E13 + EO13I  -E23 + EO23I  scal + -EOI          EO3       -EO123          EO1            0          EO2            0            0  E123 + EO123I  -E1 + -EO1I         EO13  -E2 + -EO2I         EO23          -EO        -EO12            0            0            0  -E12 + EO12I       -EO123          EO1          EO2            0         EO12 ;
         EOI          -EI        -EO1I        -EO2I        -EO3I           EO         -E1I         -E2I         -E3I         scal       -EO12I       -EO13I          EO1       -EO23I          EO2          EO3         E12I         E13I          -E1         E23I          -E2          -E3       EO123I        -EO12        -EO13        -EO23        E123I         -E12         -E13         -E23       -EO123         E123 ;
         E12        -E12I          -E2           E1         E123        -EO12          E2I         -E1I       -E123I        EO12I         scal          E23          EO2         -E13         -EO1       -EO123          -EI        -E23I        -EO2I         E13I         EO1I       EO123I           E3          -EO        -EO23         EO13         -E3I          EOI        EO23I       -EO13I         -EO3         EO3I ;
         E13        -E13I          -E3        -E123           E1        -EO13          E3I        E123I         -E1I        EO13I         -E23         scal          EO3          E12        EO123         -EO1         E23I          -EI        -EO3I        -E12I      -EO123I         EO1I          -E2         EO23          -EO        -EO12          E2I       -EO23I          EOI        EO12I          EO2        -EO2I ;
         E1I            0          -EI        -E12I        -E13I   E1 + -EO1I            0            0            0         -E1I         -E2I         -E3I   scal + EOI       -E123I  E12 + EO12I  E13 + EO13I            0            0          -EI            0        -E12I        -E13I         E23I   -E2 + EO2I   -E3 + EO3I  -E123 + EO123I            0          E2I          E3I        E123I  -E23 + -EO23I         E23I ;
         E23        -E23I         E123          -E3           E2        -EO23       -E123I          E3I         -E2I        EO23I          E13         -E12       -EO123         scal          EO3         -EO2        -E13I         E12I       EO123I          -EI        -EO3I         EO2I           E1        -EO13         EO12          -EO         -E1I        EO13I       -EO12I          EOI         -EO1         EO1I ;
         E2I            0         E12I          -EI        -E23I   E2 + -EO2I            0            0            0         -E2I          E1I        E123I  -E12 + -EO12I         -E3I   scal + EOI  E23 + EO23I            0            0         E12I            0          -EI        -E23I        -E13I   E1 + -EO1I  E123 + -EO123I   -E3 + EO3I            0         -E1I       -E123I          E3I  E13 + EO13I        -E13I ;
         E3I            0         E13I         E23I          -EI   E3 + -EO3I            0            0            0         -E3I       -E123I          E1I  -E13 + -EO13I          E2I  -E23 + -EO23I   scal + EOI            0            0         E13I            0         E23I          -EI         E12I  -E123 + EO123I   E1 + -EO1I   E2 + -EO2I            0        E123I         -E1I         -E2I  -E12 + -EO12I         E12I ;
        EO12  E12 + -EO12I         -EO2          EO1        EO123            0    E2 + EO2I  -E1 + -EO1I  -E123 + -EO123I         EO12           EO         EO23            0        -EO13            0            0  scal + -EOI  E23 + -EO23I          EO2  -E13 + EO13I         -EO1       -EO123          EO3            0            0            0  -E3 + -EO3I           EO         EO23        -EO13            0         -EO3 ;
        EO13  E13 + -EO13I         -EO3       -EO123          EO1            0    E3 + EO3I  E123 + EO123I  -E1 + -EO1I         EO13        -EO23           EO            0         EO12            0            0  -E23 + EO23I  scal + -EOI          EO3  E12 + -EO12I        EO123         -EO1         -EO2            0            0            0    E2 + EO2I        -EO23           EO         EO12            0          EO2 ;
        EO1I          E1I         -EOI       -EO12I       -EO13I          EO1           EI         E12I         E13I          -E1        -EO2I        -EO3I           EO      -EO123I         EO12         EO13         -E2I         -E3I         scal       -E123I          E12          E13        EO23I         -EO2         -EO3       -EO123        -E23I           E2           E3         E123        -EO23         -E23 ;
        EO23  E23 + -EO23I        EO123         -EO3          EO2            0  -E123 + -EO123I    E3 + EO3I  -E2 + -EO2I         EO23         EO13        -EO12            0           EO            0            0  E13 + -EO13I  -E12 + EO12I       -EO123  scal + -EOI          EO3         -EO2          EO1            0            0            0  -E1 + -EO1I         EO13        -EO12           EO            0         -EO1 ;
        EO2I          E2I        EO12I         -EOI       -EO23I          EO2        -E12I           EI         E23I          -E2         EO1I       EO123I        -EO12        -EO3I           EO         EO23          E1I        E123I         -E12         -E3I         scal          E23       -EO13I          EO1        EO123         -EO3         E13I          -E1        -E123           E3         EO13          E13 ;
        EO3I          E3I        EO13I        EO23I         -EOI          EO3        -E13I        -E23I           EI          -E3      -EO123I         EO1I        -EO13         EO2I        -EO23           EO       -E123I          E1I         -E13          E2I         -E23         scal        EO12I       -EO123          EO1          EO2        -E12I         E123          -E1          -E2        -EO12         -E12 ;
        E123       -E123I          E23         -E13          E12        EO123        -E23I         E13I        -E12I      -EO123I           E3          -E2         EO23           E1        -EO13         EO12         -E3I          E2I       -EO23I         -E1I        EO13I       -EO12I         scal          EO3         -EO2          EO1          -EI        -EO3I         EO2I        -EO1I           EO         -EOI ;
        E12I            0          E2I         -E1I       -E123I  E12 + EO12I            0            0            0        -E12I           EI         E23I   -E2 + EO2I        -E13I   E1 + -EO1I  E123 + -EO123I            0            0          E2I            0         -E1I       -E123I         -E3I   scal + EOI  E23 + EO23I  -E13 + -EO13I            0          -EI        -E23I         E13I   E3 + -EO3I         -E3I ;
        E13I            0          E3I        E123I         -E1I  E13 + EO13I            0            0            0        -E13I        -E23I           EI   -E3 + EO3I         E12I  -E123 + EO123I   E1 + -EO1I            0            0          E3I            0        E123I         -E1I          E2I  -E23 + -EO23I   scal + EOI  E12 + EO12I            0         E23I          -EI        -E12I   -E2 + EO2I          E2I ;
        E23I            0       -E123I          E3I         -E2I  E23 + EO23I            0            0            0        -E23I         E13I        -E12I  E123 + -EO123I           EI   -E3 + EO3I   E2 + -EO2I            0            0       -E123I            0          E3I         -E2I         -E1I  E13 + EO13I  -E12 + -EO12I   scal + EOI            0        -E13I         E12I          -EI   E1 + -EO1I         -E1I ;
       EO123  -E123 + -EO123I         EO23        -EO13         EO12            0  E23 + -EO23I  -E13 + EO13I  E12 + -EO12I        EO123          EO3         -EO2            0          EO1            0            0  -E3 + -EO3I    E2 + EO2I        -EO23  -E1 + -EO1I         EO13        -EO12           EO            0            0            0  scal + -EOI          EO3         -EO2          EO1            0          -EO ;
       EO12I        -E12I         EO2I        -EO1I      -EO123I         EO12          E2I         -E1I       -E123I          E12          EOI        EO23I         -EO2       -EO13I          EO1        EO123          -EI        -E23I           E2         E13I          -E1        -E123        -EO3I           EO         EO23        -EO13         -E3I         scal          E23         -E13          EO3          -E3 ;
       EO13I        -E13I         EO3I       EO123I        -EO1I         EO13          E3I        E123I         -E1I          E13       -EO23I          EOI         -EO3        EO12I       -EO123          EO1         E23I          -EI           E3        -E12I         E123          -E1         EO2I        -EO23           EO         EO12          E2I         -E23         scal          E12         -EO2           E2 ;
       EO23I        -E23I      -EO123I         EO3I        -EO2I         EO23       -E123I          E3I         -E2I          E23        EO13I       -EO12I        EO123          EOI         -EO3          EO2        -E13I         E12I        -E123          -EI           E3          -E2        -EO1I         EO13        -EO12           EO         -E1I          E13         -E12         scal          EO1          -E1 ;
       E123I            0        -E23I         E13I        -E12I  E123 + -EO123I            0            0            0       -E123I          E3I         -E2I  E23 + EO23I          E1I  -E13 + -EO13I  E12 + EO12I            0            0        -E23I            0         E13I        -E12I          -EI   E3 + -EO3I   -E2 + EO2I   E1 + -EO1I            0         -E3I          E2I         -E1I   scal + EOI          -EI ;
      EO123I        E123I       -EO23I        EO13I       -EO12I        EO123         E23I        -E13I         E12I        -E123         EO3I        -EO2I         EO23         EO1I        -EO13         EO12          E3I         -E2I          E23          E1I         -E13          E12         -EOI          EO3         -EO2          EO1           EI          -E3           E2          -E1           EO         scal ;

      ];

            R = CGA(Am*B.m);
        end

        function [scal,  EO, E1, E2, E3, EI,  EO1, EO2, EO3, EOI, E12, E13, E1I, E23, E2I, E3I,  EO12, EO13, EO1I, EO23, EO2I, EO3I, E123, E12I, E13I, E23I, EO123, EO12I, EO13I, EO23I, E123I, EO123I] = decompose_(A)
            scal = A.m(1);

            EO = A.m(2); E1 = A.m(3); E2 = A.m(4); E3 = A.m(5); EI = A.m(6);

            EO1 = A.m(7); EO2 = A.m(8); EO3 = A.m(9); EOI = A.m(10);
            E12 = A.m(11); E13 = A.m(12); E1I = A.m(13);
            E23 = A.m(14); E2I = A.m(15);
	        E3I = A.m(16);
	    
            EO12 = A.m(17); EO13 = A.m(18); EO1I = A.m(19);
            EO23 = A.m(20); EO2I = A.m(21);
            EO3I = A.m(22);
            E123 = A.m(23); E12I = A.m(24); 
            E13I = A.m(25);
            E23I = A.m(26);

            EO123 = A.m(27); EO12I = A.m(28);
            EO13I = A.m(29);
            EO23I = A.m(30);
            E123I = A.m(31);
            
            EO123I = A.m(32);
        end

        % ***** The Outer Product *****

        
        function R = outer_(A, B)

             [scal,  EO, E1, E2, E3, EI,  EO1, EO2, EO3, EOI, E12, E13, E1I, E23, E2I, E3I,  EO12, EO13, EO1I, EO23, EO2I, EO3I, E123, E12I, E13I, E23I, EO123, EO12I, EO13I, EO23I, E123I, EO123I] = decompose_(A);
rm = [
   scal      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0 ;
    EO   scal      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0 ;
    E1      0   scal      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0 ;
    E2      0      0   scal      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0 ;
    E3      0      0      0   scal      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0 ;
    EI      0      0      0      0   scal      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0 ;
   EO1    -E1     EO      0      0      0   scal      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0 ;
   EO2    -E2      0     EO      0      0      0   scal      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0 ;
   EO3    -E3      0      0     EO      0      0      0   scal      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0 ;
   EOI    -EI      0      0      0     EO      0      0      0   scal      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0 ;
   E12      0    -E2     E1      0      0      0      0      0      0   scal      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0 ;
   E13      0    -E3      0     E1      0      0      0      0      0      0   scal      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0 ;
   E1I      0    -EI      0      0     E1      0      0      0      0      0      0   scal      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0 ;
   E23      0      0    -E3     E2      0      0      0      0      0      0      0      0   scal      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0 ;
   E2I      0      0    -EI      0     E2      0      0      0      0      0      0      0      0   scal      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0 ;
   E3I      0      0      0    -EI     E3      0      0      0      0      0      0      0      0      0   scal      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0 ;
  EO12    E12   -EO2    EO1      0      0     E2    -E1      0      0     EO      0      0      0      0      0   scal      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0 ;
  EO13    E13   -EO3      0    EO1      0     E3      0    -E1      0      0     EO      0      0      0      0      0   scal      0      0      0      0      0      0      0      0      0      0      0      0      0      0 ;
  EO1I    E1I   -EOI      0      0    EO1     EI      0      0    -E1      0      0     EO      0      0      0      0      0   scal      0      0      0      0      0      0      0      0      0      0      0      0      0 ;
  EO23    E23      0   -EO3    EO2      0      0     E3    -E2      0      0      0      0     EO      0      0      0      0      0   scal      0      0      0      0      0      0      0      0      0      0      0      0 ;
  EO2I    E2I      0   -EOI      0    EO2      0     EI      0    -E2      0      0      0      0     EO      0      0      0      0      0   scal      0      0      0      0      0      0      0      0      0      0      0 ;
  EO3I    E3I      0      0   -EOI    EO3      0      0     EI    -E3      0      0      0      0      0     EO      0      0      0      0      0   scal      0      0      0      0      0      0      0      0      0      0 ;
  E123      0    E23   -E13    E12      0      0      0      0      0     E3    -E2      0     E1      0      0      0      0      0      0      0      0   scal      0      0      0      0      0      0      0      0      0 ;
  E12I      0    E2I   -E1I      0    E12      0      0      0      0     EI      0    -E2      0     E1      0      0      0      0      0      0      0      0   scal      0      0      0      0      0      0      0      0 ;
  E13I      0    E3I      0   -E1I    E13      0      0      0      0      0     EI    -E3      0      0     E1      0      0      0      0      0      0      0      0   scal      0      0      0      0      0      0      0 ;
  E23I      0      0    E3I   -E2I    E23      0      0      0      0      0      0      0     EI    -E3     E2      0      0      0      0      0      0      0      0      0   scal      0      0      0      0      0      0 ;
 EO123  -E123   EO23  -EO13   EO12      0    E23   -E13    E12      0    EO3   -EO2      0    EO1      0      0    -E3     E2      0    -E1      0      0     EO      0      0      0   scal      0      0      0      0      0 ;
 EO12I  -E12I   EO2I  -EO1I      0   EO12    E2I   -E1I      0    E12    EOI      0   -EO2      0    EO1      0    -EI      0     E2      0    -E1      0      0     EO      0      0      0   scal      0      0      0      0 ;
 EO13I  -E13I   EO3I      0  -EO1I   EO13    E3I      0   -E1I    E13      0    EOI   -EO3      0      0    EO1      0    -EI     E3      0      0    -E1      0      0     EO      0      0      0   scal      0      0      0 ;
 EO23I  -E23I      0   EO3I  -EO2I   EO23      0    E3I   -E2I    E23      0      0      0    EOI   -EO3    EO2      0      0      0    -EI     E3    -E2      0      0      0     EO      0      0      0   scal      0      0 ;
 E123I      0  -E23I   E13I  -E12I   E123      0      0      0      0    E3I   -E2I    E23    E1I   -E13    E12      0      0      0      0      0      0    -EI     E3    -E2     E1      0      0      0      0   scal      0 ;
 EO123I  E123I  -EO23I  EO13I  -EO12I  EO123   E23I  -E13I   E12I  -E123   EO3I  -EO2I   EO23   EO1I  -EO13   EO12    E3I   -E2I    E23    E1I   -E13    E12   -EOI    EO3   -EO2    EO1     EI    -E3     E2    -E1     EO   scal ;

 ];

                R = CGA(rm*B.m);
        end

        function R = divide_(A, B)
            R = A * inverse_(B);
        end

        function R = leftcontraction_(A, B)

            [scal,  EO, E1, E2, E3, EI,  EO1, EO2, EO3, EOI, E12, E13, E1I, E23, E2I, E3I,  EO12, EO13, EO1I, EO23, EO2I, EO3I, E123, E12I, E13I, E23I, EO123, EO12I, EO13I, EO23I, E123I, EO123I] = decompose_(A);
            
            M = [
 scal   -EI      E1      E2      E3     -EO     -E1I    -E2I    -E3I     EOI    -E12    -E13    -EO1    -E23    -EO2    -EO3     E12I    E13I    EO1I    E23I    EO2I    EO3I   -E123    EO12    EO13    EO23    E123I  -EO12I  -EO13I  -EO23I   EO123  -EO123I ;
0       scal   0      0      0      0      -E1     -E2     -E3      EO     0      0      0      0      0      0      -E12    -E13    -EO1    -E23    -EO2    -EO3    0      0      0      0       E123   -EO12   -EO13   -EO23   0       EO123  ;
0      0       scal   0      0      0      -EI     0      0      0      -E2     -E3      EO     0      0      0       E2I     E3I    -EOI    0      0      0      -E23    -EO2    -EO3    0       E23I    EO2I    EO3I   0      -EO23    EO23I  ;
0      0      0       scal   0      0      0      -EI     0      0       E1     0      0      -E3      EO     0      -E1I    0      0       E3I    -EOI    0       E13     EO1    0      -EO3    -E13I   -EO1I   0       EO3I    EO13   -EO13I  ;
0      0      0      0       scal   0      0      0      -EI     0      0       E1     0       E2     0       EO     0      -E1I    0      -E2I    0      -EOI    -E12    0       EO1     EO2     E12I   0      -EO1I   -EO2I   -EO12    EO12I  ;
0      0      0      0      0       scal   0      0      0      -EI     0      0       E1     0       E2      E3     0      0      -E1I    0      -E2I    -E3I    0      -E12    -E13    -E23    0       E12I    E13I    E23I   -E123    E123I  ;
0      0      0      0      0      0       scal   0      0      0      0      0      0      0      0      0       E2      E3     -EO     0      0      0      0      0      0      0      -E23    -EO2    -EO3    0      0       EO23   ;
0      0      0      0      0      0      0       scal   0      0      0      0      0      0      0      0      -E1     0      0       E3     -EO     0      0      0      0      0       E13     EO1    0      -EO3    0      -EO13   ;
0      0      0      0      0      0      0      0       scal   0      0      0      0      0      0      0      0      -E1     0      -E2     0      -EO     0      0      0      0      -E12    0       EO1     EO2    0       EO12   ;
0      0      0      0      0      0      0      0      0       scal   0      0      0      0      0      0      0      0      -E1     0      -E2     -E3     0      0      0      0      0      -E12    -E13    -E23    0       E123   ;
0      0      0      0      0      0      0      0      0      0       scal   0      0      0      0      0      -EI     0      0      0      0      0       E3     -EO     0      0      -E3I     EOI    0      0      -EO3     EO3I   ;
0      0      0      0      0      0      0      0      0      0      0       scal   0      0      0      0      0      -EI     0      0      0      0      -E2     0      -EO     0       E2I    0       EOI    0       EO2    -EO2I   ;
0      0      0      0      0      0      0      0      0      0      0      0       scal   0      0      0      0      0      -EI     0      0      0      0      -E2     -E3     0      0       E2I     E3I    0      -E23     E23I   ;
0      0      0      0      0      0      0      0      0      0      0      0      0       scal   0      0      0      0      0      -EI     0      0       E1     0      0      -EO     -E1I    0      0       EOI    -EO1     EO1I   ;
0      0      0      0      0      0      0      0      0      0      0      0      0      0       scal   0      0      0      0      0      -EI     0      0       E1     0      -E3     0      -E1I    0       E3I     E13    -E13I   ;
0      0      0      0      0      0      0      0      0      0      0      0      0      0      0       scal   0      0      0      0      0      -EI     0      0       E1      E2     0      0      -E1I    -E2I    -E12     E12I   ;
0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0       scal   0      0      0      0      0      0      0      0      0      -E3      EO     0      0      0      -EO3    ;
0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0       scal   0      0      0      0      0      0      0      0       E2     0       EO     0      0       EO2    ;
0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0       scal   0      0      0      0      0      0      0      0       E2      E3     0      0      -E23    ;
0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0       scal   0      0      0      0      0      0      -E1     0      0       EO     0      -EO1    ;
0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0       scal   0      0      0      0      0      0      -E1     0       E3     0       E13    ;
0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0       scal   0      0      0      0      0      0      -E1     -E2     0      -E12    ;
0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0       scal   0      0      0      -EI     0      0      0       EO     -EOI    ;
0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0       scal   0      0      0      -EI     0      0       E3     -E3I    ;
0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0       scal   0      0      0      -EI     0      -E2      E2I    ;
0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0       scal   0      0      0      -EI      E1     -E1I    ;
0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0       scal   0      0      0      0      -EO     ;
0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0       scal   0      0      0      -E3     ;
0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0       scal   0      0       E2     ;
0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0       scal   0      -E1     ;
0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0       scal   -EI     ;
0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0       scal   ;
		  ];
            R = CGA(M*B.m);
        end

        function R = rightcontraction_(A, B)
            [scal,  EO, E1, E2, E3, EI,  EO1, EO2, EO3, EOI, E12, E13, E1I, E23, E2I, E3I,  EO12, EO13, EO1I, EO23, EO2I, EO3I, E123, E12I, E13I, E23I, EO123, EO12I, EO13I, EO23I, E123I, EO123I] = decompose_(A);

	    M = [
 scal   -EI      E1      E2      E3     -EO     -E1I    -E2I    -E3I     EOI    -E12    -E13    -EO1    -E23    -EO2    -EO3     E12I    E13I    EO1I    E23I    EO2I    EO3I   -E123    EO12    EO13    EO23    E123I  -EO12I  -EO13I  -EO23I   EO123  -EO123I ;
 EO     -EOI     EO1     EO2     EO3    0      -EO1I   -EO2I   -EO3I   0      -EO12   -EO13   0      -EO23   0      0       EO12I   EO13I  0       EO23I  0      0      -EO123  0      0      0       EO123I 0      0      0      0      0      ;
 E1     -E1I    0       E12     E13     EO1    0      -E12I   -E13I   -EO1I   0      0      0      -E123    EO12    EO13   0      0      0       E123I  -EO12I  -EO13I  0      0      0      -EO123  0      0      0       EO123I 0      0      ;
 E2     -E2I    -E12    0       E23     EO2     E12I   0      -E23I   -EO2I   0       E123   -EO12   0      0       EO23   0      -E123I   EO12I  0      0      -EO23I  0      0       EO123  0      0      0      -EO123I 0      0      0      ;
 E3     -E3I    -E13    -E23    0       EO3     E13I    E23I   0      -EO3I   -E123   0      -EO13   0      -EO23   0       E123I  0       EO13I  0       EO23I  0      0      -EO123  0      0      0       EO123I 0      0      0      0      ;
 EI     0      -E1I    -E2I    -E3I     EOI    0      0      0      0      -E12I   -E13I   -EO1I   -E23I   -EO2I   -EO3I   0      0      0      0      0      0       E123I  -EO12I  -EO13I  -EO23I  0      0      0      0       EO123I 0      ;
 EO1    -EO1I   0       EO12    EO13   0      0      -EO12I  -EO13I  0      0      0      0      -EO123  0      0      0      0      0       EO123I 0      0      0      0      0      0      0      0      0      0      0      0      ;
 EO2    -EO2I   -EO12   0       EO23   0       EO12I  0      -EO23I  0      0       EO123  0      0      0      0      0      -EO123I 0      0      0      0      0      0      0      0      0      0      0      0      0      0      ;
 EO3    -EO3I   -EO13   -EO23   0      0       EO13I   EO23I  0      0      -EO123  0      0      0      0      0       EO123I 0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      ;
 EOI    0      -EO1I   -EO2I   -EO3I   0      0      0      0      0      -EO12I  -EO13I  0      -EO23I  0      0      0      0      0      0      0      0       EO123I 0      0      0      0      0      0      0      0      0      ;
 E12    -E12I   0      0       E123   -EO12   0      0      -E123I   EO12I  0      0      0      0      0      -EO123  0      0      0      0      0       EO123I 0      0      0      0      0      0      0      0      0      0      ;
 E13    -E13I   0      -E123   0      -EO13   0       E123I  0       EO13I  0      0      0      0       EO123  0      0      0      0      0      -EO123I 0      0      0      0      0      0      0      0      0      0      0      ;
 E1I    0      0      -E12I   -E13I   -EO1I   0      0      0      0      0      0      0      -E123I   EO12I   EO13I  0      0      0      0      0      0      0      0      0       EO123I 0      0      0      0      0      0      ;
 E23    -E23I    E123   0      0      -EO23   -E123I  0      0       EO23I  0      0      -EO123  0      0      0      0      0       EO123I 0      0      0      0      0      0      0      0      0      0      0      0      0      ;
 E2I    0       E12I   0      -E23I   -EO2I   0      0      0      0      0       E123I  -EO12I  0      0       EO23I  0      0      0      0      0      0      0      0      -EO123I 0      0      0      0      0      0      0      ;
 E3I    0       E13I    E23I   0      -EO3I   0      0      0      0      -E123I  0      -EO13I  0      -EO23I  0      0      0      0      0      0      0      0       EO123I 0      0      0      0      0      0      0      0      ;
 EO12   -EO12I  0      0       EO123  0      0      0      -EO123I 0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      ;
 EO13   -EO13I  0      -EO123  0      0      0       EO123I 0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      ;
 EO1I   0      0      -EO12I  -EO13I  0      0      0      0      0      0      0      0      -EO123I 0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      ;
 EO23   -EO23I   EO123  0      0      0      -EO123I 0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      ;
 EO2I   0       EO12I  0      -EO23I  0      0      0      0      0      0       EO123I 0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      ;
 EO3I   0       EO13I   EO23I  0      0      0      0      0      0      -EO123I 0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      ;
 E123   -E123I  0      0      0       EO123  0      0      0      -EO123I 0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      ;
 E12I   0      0      0      -E123I   EO12I  0      0      0      0      0      0      0      0      0      -EO123I 0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      ;
 E13I   0      0       E123I  0       EO13I  0      0      0      0      0      0      0      0       EO123I 0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      ;
 E23I   0      -E123I  0      0       EO23I  0      0      0      0      0      0      -EO123I 0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      ;
 EO123  -EO123I 0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      ;
 EO12I  0      0      0      -EO123I 0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      ;
 EO13I  0      0       EO123I 0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      ;
 EO23I  0      -EO123I 0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      ;
 E123I  0      0      0      0      -EO123I 0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      ;
 EO123I 0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      ;

	    ];
            R = CGA(M*B.m);
        end

        function R = inner_(A, B)
 	    [scal,  EO, E1, E2, E3, EI,  EO1, EO2, EO3, EOI, E12, E13, E1I, E23, E2I, E3I,  EO12, EO13, EO1I, EO23, EO2I, EO3I, E123, E12I, E13I, E23I, EO123, EO12I, EO13I, EO23I, E123I, EO123I] = decompose_(A);

	       M = [
        scal          -EI           E1           E2           E3          -EO         -E1I         -E2I         -E3I          EOI         -E12         -E13         -EO1         -E23         -EO2         -EO3         E12I         E13I         EO1I         E23I         EO2I         EO3I        -E123         EO12         EO13         EO23        E123I       -EO12I       -EO13I       -EO23I        EO123      -EO123I ;
          EO  scal + -EOI          EO1          EO2          EO3            0  -E1 + -EO1I  -E2 + -EO2I  -E3 + -EO3I           EO        -EO12        -EO13            0        -EO23            0            0  -E12 + EO12I  -E13 + EO13I         -EO1  -E23 + EO23I         -EO2         -EO3       -EO123            0            0            0  E123 + EO123I        -EO12        -EO13        -EO23            0        EO123 ;
          E1         -E1I         scal          E12          E13          EO1          -EI        -E12I        -E13I        -EO1I          -E2          -E3           EO        -E123         EO12         EO13          E2I          E3I         -EOI        E123I       -EO12I       -EO13I         -E23         -EO2         -EO3       -EO123         E23I         EO2I         EO3I       EO123I        -EO23        EO23I ;
          E2         -E2I         -E12         scal          E23          EO2         E12I          -EI        -E23I        -EO2I           E1         E123        -EO12          -E3           EO         EO23         -E1I       -E123I        EO12I          E3I         -EOI       -EO23I          E13          EO1        EO123         -EO3        -E13I        -EO1I      -EO123I         EO3I         EO13       -EO13I ;
          E3         -E3I         -E13         -E23         scal          EO3         E13I         E23I          -EI        -EO3I        -E123           E1        -EO13           E2        -EO23           EO        E123I         -E1I        EO13I         -E2I        EO23I         -EOI         -E12       -EO123          EO1          EO2         E12I       EO123I        -EO1I        -EO2I        -EO12        EO12I ;
          EI            0         -E1I         -E2I         -E3I   scal + EOI            0            0            0          -EI        -E12I        -E13I   E1 + -EO1I        -E23I   E2 + -EO2I   E3 + -EO3I            0            0         -E1I            0         -E2I         -E3I        E123I  -E12 + -EO12I  -E13 + -EO13I  -E23 + -EO23I            0         E12I         E13I         E23I  -E123 + EO123I        E123I ;
         EO1        -EO1I            0         EO12         EO13            0  scal + -EOI  -E12 + -EO12I  -E13 + -EO13I         -EO1         -EO2         -EO3            0       -EO123            0            0    E2 + EO2I    E3 + EO3I          -EO  -E123 + EO123I         EO12         EO13        -EO23            0            0            0  -E23 + EO23I         -EO2         -EO3        EO123            0         EO23 ;
         EO2        -EO2I        -EO12            0         EO23            0  E12 + EO12I  scal + -EOI  -E23 + -EO23I         -EO2          EO1        EO123            0         -EO3            0            0  -E1 + -EO1I  E123 + -EO123I        -EO12    E3 + EO3I          -EO         EO23         EO13            0            0            0  E13 + -EO13I          EO1       -EO123         -EO3            0        -EO13 ;
         EO3        -EO3I        -EO13        -EO23            0            0  E13 + EO13I  E23 + EO23I  scal + -EOI         -EO3       -EO123          EO1            0          EO2            0            0  -E123 + EO123I  -E1 + -EO1I        -EO13  -E2 + -EO2I        -EO23          -EO        -EO12            0            0            0  -E12 + EO12I        EO123          EO1          EO2            0         EO12 ;
         EOI            0        -EO1I        -EO2I        -EO3I            0            0            0            0         scal       -EO12I       -EO13I            0       -EO23I            0            0            0            0          -E1            0          -E2          -E3       EO123I            0            0            0            0         -E12         -E13         -E23            0         E123 ;
         E12        -E12I            0            0         E123        -EO12         -E2I         -E1I       -E123I        EO12I         scal         -E23         -EO2         -E13         -EO1       -EO123          -EI         E23I         EO2I         E13I         EO1I       EO123I           E3          -EO         EO23         EO13         -E3I          EOI       -EO23I       -EO13I         -EO3         EO3I ;
         E13        -E13I            0        -E123            0        -EO13         -E3I        E123I         -E1I        EO13I          E23         scal         -EO3          E12        EO123         -EO1        -E23I          -EI         EO3I        -E12I      -EO123I         EO1I          -E2        -EO23          -EO        -EO12          E2I        EO23I          EOI        EO12I          EO2        -EO2I ;
         E1I            0            0        -E12I        -E13I        -EO1I            0            0            0         -E1I          E2I          E3I  scal + -EOI       -E123I  E12 + EO12I  E13 + EO13I            0            0          -EI            0        -E12I        -E13I        -E23I  -E2 + -EO2I  -E3 + -EO3I  -E123 + EO123I            0          E2I          E3I        E123I  -E23 + EO23I         E23I ;
         E23        -E23I         E123            0            0        -EO23       -E123I         -E3I         -E2I        EO23I         -E13         -E12       -EO123         scal         -EO3         -EO2         E13I         E12I       EO123I          -EI         EO3I         EO2I           E1         EO13         EO12          -EO         -E1I       -EO13I       -EO12I          EOI         -EO1         EO1I ;
         E2I            0         E12I            0        -E23I        -EO2I            0            0            0         -E2I         -E1I        E123I  -E12 + -EO12I          E3I  scal + -EOI  E23 + EO23I            0            0         E12I            0          -EI        -E23I         E13I    E1 + EO1I  E123 + -EO123I  -E3 + -EO3I            0         -E1I       -E123I          E3I  E13 + -EO13I        -E13I ;
         E3I            0         E13I         E23I            0        -EO3I            0            0            0         -E3I       -E123I         -E1I  -E13 + -EO13I         -E2I  -E23 + -EO23I  scal + -EOI            0            0         E13I            0         E23I          -EI        -E12I  -E123 + EO123I    E1 + EO1I    E2 + EO2I            0        E123I         -E1I         -E2I  -E12 + EO12I         E12I ;
        EO12       -EO12I            0            0        EO123            0        -EO2I        -EO1I  -E123 + -EO123I         EO12            0        -EO23            0        -EO13            0            0  scal + -EOI  E23 + EO23I          EO2  E13 + EO13I          EO1       -EO123          EO3            0            0            0  -E3 + -EO3I           EO         EO23         EO13            0         -EO3 ;
        EO13       -EO13I            0       -EO123            0            0        -EO3I  E123 + EO123I        -EO1I         EO13         EO23            0            0         EO12            0            0  -E23 + -EO23I  scal + -EOI          EO3  -E12 + -EO12I        EO123          EO1         -EO2            0            0            0    E2 + EO2I        -EO23           EO        -EO12            0          EO2 ;
        EO1I            0            0       -EO12I       -EO13I            0            0            0            0            0         EO2I         EO3I            0      -EO123I            0            0            0            0         scal            0         -E12         -E13       -EO23I            0            0            0            0           E2           E3        -E123            0         -E23 ;
        EO23       -EO23I        EO123            0            0            0  -E123 + -EO123I        -EO3I        -EO2I         EO23        -EO13        -EO12            0            0            0            0  E13 + EO13I  E12 + EO12I       -EO123  scal + -EOI          EO3          EO2          EO1            0            0            0  -E1 + -EO1I         EO13         EO12           EO            0         -EO1 ;
        EO2I            0        EO12I            0       -EO23I            0            0            0            0            0        -EO1I       EO123I            0         EO3I            0            0            0            0          E12            0         scal         -E23        EO13I            0            0            0            0          -E1         E123           E3            0          E13 ;
        EO3I            0        EO13I        EO23I            0            0            0            0            0            0      -EO123I        -EO1I            0        -EO2I            0            0            0            0          E13            0          E23         scal       -EO12I            0            0            0            0        -E123          -E1          -E2            0         -E12 ;
        E123       -E123I            0            0            0        EO123        -E23I        -E13I        -E12I      -EO123I            0            0         EO23            0         EO13         EO12         -E3I         -E2I       -EO23I         -E1I       -EO13I       -EO12I         scal          EO3          EO2          EO1          -EI        -EO3I        -EO2I        -EO1I           EO         -EOI ;
        E12I            0            0            0       -E123I        EO12I            0            0            0        -E12I            0         E23I         EO2I         E13I         EO1I  E123 + -EO123I            0            0         -E2I            0         -E1I       -E123I         -E3I   scal + EOI  -E23 + EO23I  -E13 + EO13I            0          -EI         E23I         E13I   E3 + -EO3I         -E3I ;
        E13I            0            0        E123I            0        EO13I            0            0            0        -E13I        -E23I            0         EO3I        -E12I  -E123 + EO123I         EO1I            0            0         -E3I            0        E123I         -E1I          E2I  E23 + -EO23I   scal + EOI  E12 + -EO12I            0        -E23I          -EI        -E12I   -E2 + EO2I          E2I ;
        E23I            0       -E123I            0            0        EO23I            0            0            0        -E23I         E13I         E12I  E123 + -EO123I            0         EO3I         EO2I            0            0       -E123I            0         -E3I         -E2I         -E1I  -E13 + EO13I  -E12 + EO12I   scal + EOI            0         E13I         E12I          -EI   E1 + -EO1I         -E1I ;
       EO123      -EO123I            0            0            0            0       -EO23I       -EO13I       -EO12I       -EO123            0            0            0            0            0            0        -EO3I        -EO2I        -EO23        -EO1I        -EO13        -EO12            0            0            0            0  scal + -EOI         -EO3         -EO2         -EO1            0          -EO ;
       EO12I            0            0            0      -EO123I            0            0            0            0            0            0        EO23I            0        EO13I            0            0            0            0            0            0            0        -E123        -EO3I            0            0            0            0         scal          E23          E13            0          -E3 ;
       EO13I            0            0       EO123I            0            0            0            0            0            0       -EO23I            0            0       -EO12I            0            0            0            0            0            0         E123            0         EO2I            0            0            0            0         -E23         scal         -E12            0           E2 ;
       EO23I            0      -EO123I            0            0            0            0            0            0            0        EO13I        EO12I            0            0            0            0            0            0        -E123            0            0            0        -EO1I            0            0            0            0          E13          E12         scal            0          -E1 ;
       E123I            0            0            0            0      -EO123I            0            0            0       -E123I            0            0       -EO23I            0       -EO13I       -EO12I            0            0        -E23I            0        -E13I        -E12I            0        -EO3I        -EO2I        -EO1I            0         -E3I         -E2I         -E1I  scal + -EOI          -EI ;
      EO123I            0            0            0            0            0            0            0            0            0            0            0            0            0            0            0            0            0            0            0            0            0            0            0            0            0            0            0            0            0            0         scal ;

	       ];
M = [
 scal   -EI      E1      E2      E3     -EO     -E1I    -E2I    -E3I     EOI    -E12    -E13    -EO1    -E23    -EO2    -EO3     E12I    E13I    EO1I    E23I    EO2I    EO3I   -E123    EO12    EO13    EO23    E123I  -EO12I  -EO13I  -EO23I   EO123  -EO123I ;
 EO     -EOI     EO1     EO2     EO3    0      -EO1I   -EO2I   -EO3I    EO     -EO12   -EO13   0      -EO23   0      0       EO12I   EO13I  -EO1     EO23I  -EO2    -EO3    -EO123  0      0      0       EO123I -EO12   -EO13   -EO23   0       EO123  ;
 E1     -E1I     scal    E12     E13     EO1    -EI     -E12I   -E13I   -EO1I   -E2     -E3      EO     -E123    EO12    EO13    E2I     E3I    -EOI     E123I  -EO12I  -EO13I  -E23    -EO2    -EO3    -EO123   E23I    EO2I    EO3I    EO123I -EO23    EO23I  ;
 E2     -E2I    -E12     scal    E23     EO2     E12I   -EI     -E23I   -EO2I    E1      E123   -EO12   -E3      EO      EO23   -E1I    -E123I   EO12I   E3I    -EOI    -EO23I   E13     EO1     EO123  -EO3    -E13I   -EO1I   -EO123I  EO3I    EO13   -EO13I  ;
 E3     -E3I    -E13    -E23     scal    EO3     E13I    E23I   -EI     -EO3I   -E123    E1     -EO13    E2     -EO23    EO      E123I  -E1I     EO13I  -E2I     EO23I  -EOI    -E12    -EO123   EO1     EO2     E12I    EO123I -EO1I   -EO2I   -EO12    EO12I  ;
 EI     0      -E1I    -E2I    -E3I     EOI    0      0      0      -EI     -E12I   -E13I   -EO1I   -E23I   -EO2I   -EO3I   0      0      -E1I    0      -E2I    -E3I     E123I  -EO12I  -EO13I  -EO23I  0       E12I    E13I    E23I    EO123I  E123I  ;
 EO1    -EO1I   0       EO12    EO13   0       scal   -EO12I  -EO13I  0      0      0      0      -EO123  0      0       E2      E3     -EO      EO123I 0      0      0      0      0      0      -E23    -EO2    -EO3    0      0       EO23   ;
 EO2    -EO2I   -EO12   0       EO23   0       EO12I   scal   -EO23I  0      0       EO123  0      0      0      0      -E1     -EO123I 0       E3     -EO     0      0      0      0      0       E13     EO1    0      -EO3    0      -EO13   ;
 EO3    -EO3I   -EO13   -EO23   0      0       EO13I   EO23I   scal   0      -EO123  0      0      0      0      0       EO123I -E1     0      -E2     0      -EO     0      0      0      0      -E12    0       EO1     EO2    0       EO12   ;
 EOI    0      -EO1I   -EO2I   -EO3I   0      0      0      0       scal   -EO12I  -EO13I  0      -EO23I  0      0      0      0      -E1     0      -E2     -E3      EO123I 0      0      0      0      -E12    -E13    -E23    0       E123   ;
 E12    -E12I   0      0       E123   -EO12   0      0      -E123I   EO12I   scal   0      0      0      0      -EO123  -EI     0      0      0      0       EO123I  E3     -EO     0      0      -E3I     EOI    0      0      -EO3     EO3I   ;
 E13    -E13I   0      -E123   0      -EO13   0       E123I  0       EO13I  0       scal   0      0       EO123  0      0      -EI     0      0      -EO123I 0      -E2     0      -EO     0       E2I    0       EOI    0       EO2    -EO2I   ;
 E1I    0      0      -E12I   -E13I   -EO1I   0      0      0      0      0      0       scal   -E123I   EO12I   EO13I  0      0      -EI     0      0      0      0      -E2     -E3      EO123I 0       E2I     E3I    0      -E23     E23I   ;
 E23    -E23I    E123   0      0      -EO23   -E123I  0      0       EO23I  0      0      -EO123   scal   0      0      0      0       EO123I -EI     0      0       E1     0      0      -EO     -E1I    0      0       EOI    -EO1     EO1I   ;
 E2I    0       E12I   0      -E23I   -EO2I   0      0      0      0      0       E123I  -EO12I  0       scal    EO23I  0      0      0      0      -EI     0      0       E1     -EO123I -E3     0      -E1I    0       E3I     E13    -E13I   ;
 E3I    0       E13I    E23I   0      -EO3I   0      0      0      0      -E123I  0      -EO13I  0      -EO23I   scal   0      0      0      0      0      -EI     0       EO123I  E1      E2     0      0      -E1I    -E2I    -E12     E12I   ;
 EO12   -EO12I  0      0       EO123  0      0      0      -EO123I 0      0      0      0      0      0      0       scal   0      0      0      0      0      0      0      0      0      -E3      EO     0      0      0      -EO3    ;
 EO13   -EO13I  0      -EO123  0      0      0       EO123I 0      0      0      0      0      0      0      0      0       scal   0      0      0      0      0      0      0      0       E2     0       EO     0      0       EO2    ;
 EO1I   0      0      -EO12I  -EO13I  0      0      0      0      0      0      0      0      -EO123I 0      0      0      0       scal   0      0      0      0      0      0      0      0       E2      E3     0      0      -E23    ;
 EO23   -EO23I   EO123  0      0      0      -EO123I 0      0      0      0      0      0      0      0      0      0      0      0       scal   0      0      0      0      0      0      -E1     0      0       EO     0      -EO1    ;
 EO2I   0       EO12I  0      -EO23I  0      0      0      0      0      0       EO123I 0      0      0      0      0      0      0      0       scal   0      0      0      0      0      0      -E1     0       E3     0       E13    ;
 EO3I   0       EO13I   EO23I  0      0      0      0      0      0      -EO123I 0      0      0      0      0      0      0      0      0      0       scal   0      0      0      0      0      0      -E1     -E2     0      -E12    ;
 E123   -E123I  0      0      0       EO123  0      0      0      -EO123I 0      0      0      0      0      0      0      0      0      0      0      0       scal   0      0      0      -EI     0      0      0       EO     -EOI    ;
 E12I   0      0      0      -E123I   EO12I  0      0      0      0      0      0      0      0      0      -EO123I 0      0      0      0      0      0      0       scal   0      0      0      -EI     0      0       E3     -E3I    ;
 E13I   0      0       E123I  0       EO13I  0      0      0      0      0      0      0      0       EO123I 0      0      0      0      0      0      0      0      0       scal   0      0      0      -EI     0      -E2      E2I    ;
 E23I   0      -E123I  0      0       EO23I  0      0      0      0      0      0      -EO123I 0      0      0      0      0      0      0      0      0      0      0      0       scal   0      0      0      -EI      E1     -E1I    ;
 EO123  -EO123I 0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0       scal   0      0      0      0      -EO     ;
 EO12I  0      0      0      -EO123I 0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0       scal   0      0      0      -E3     ;
 EO13I  0      0       EO123I 0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0       scal   0      0       E2     ;
 EO23I  0      -EO123I 0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0       scal   0      -E1     ;
 E123I  0      0      0      0      -EO123I 0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0       scal   -EI     ;
 EO123I 0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0       scal   ;

];
            R = CGA(M*B.m);
        end

	function R= productleftexpand_(A)
	 [scal,  EO, E1, E2, E3, EI,  EO1, EO2, EO3, EOI, E12, E13, E1I, E23, E2I, E3I,  EO12, EO13, EO1I, EO23, EO2I, EO3I, E123, E12I, E13I, E23I, EO123, EO12I, EO13I, EO23I, E123I, EO123I] = decompose_(A);

R = [
        scal          -EI           E1           E2           E3          -EO         -E1I         -E2I         -E3I          EOI         -E12         -E13         -EO1         -E23         -EO2         -EO3         E12I         E13I         EO1I         E23I         EO2I         EO3I        -E123         EO12         EO13         EO23        E123I       -EO12I       -EO13I       -EO23I        EO123      -EO123I ;
          EO  scal + -EOI          EO1          EO2          EO3            0  -E1 + -EO1I  -E2 + -EO2I  -E3 + -EO3I           EO        -EO12        -EO13            0        -EO23            0            0  -E12 + EO12I  -E13 + EO13I         -EO1  -E23 + EO23I         -EO2         -EO3       -EO123            0            0            0  E123 + EO123I        -EO12        -EO13        -EO23            0        EO123 ;
          E1         -E1I         scal          E12          E13          EO1          -EI        -E12I        -E13I        -EO1I          -E2          -E3           EO        -E123         EO12         EO13          E2I          E3I         -EOI        E123I       -EO12I       -EO13I         -E23         -EO2         -EO3       -EO123         E23I         EO2I         EO3I       EO123I        -EO23        EO23I ;
          E2         -E2I         -E12         scal          E23          EO2         E12I          -EI        -E23I        -EO2I           E1         E123        -EO12          -E3           EO         EO23         -E1I       -E123I        EO12I          E3I         -EOI       -EO23I          E13          EO1        EO123         -EO3        -E13I        -EO1I      -EO123I         EO3I         EO13       -EO13I ;
          E3         -E3I         -E13         -E23         scal          EO3         E13I         E23I          -EI        -EO3I        -E123           E1        -EO13           E2        -EO23           EO        E123I         -E1I        EO13I         -E2I        EO23I         -EOI         -E12       -EO123          EO1          EO2         E12I       EO123I        -EO1I        -EO2I        -EO12        EO12I ;
          EI            0         -E1I         -E2I         -E3I   scal + EOI            0            0            0          -EI        -E12I        -E13I   E1 + -EO1I        -E23I   E2 + -EO2I   E3 + -EO3I            0            0         -E1I            0         -E2I         -E3I        E123I  -E12 + -EO12I  -E13 + -EO13I  -E23 + -EO23I            0         E12I         E13I         E23I  -E123 + EO123I        E123I ;
         EO1  -E1 + -EO1I           EO         EO12         EO13            0  scal + -EOI  E12 + -EO12I  E13 + -EO13I          EO1         -EO2         -EO3            0       -EO123            0            0    E2 + EO2I    E3 + EO3I          -EO  E123 + EO123I        -EO12        -EO13        -EO23            0            0            0  -E23 + EO23I         -EO2         -EO3       -EO123            0         EO23 ;
         EO2  -E2 + -EO2I        -EO12           EO         EO23            0  -E12 + EO12I  scal + -EOI  E23 + -EO23I          EO2          EO1        EO123            0         -EO3            0            0  -E1 + -EO1I  -E123 + -EO123I         EO12    E3 + EO3I          -EO        -EO23         EO13            0            0            0  E13 + -EO13I          EO1        EO123         -EO3            0        -EO13 ;
         EO3  -E3 + -EO3I        -EO13        -EO23           EO            0  -E13 + EO13I  -E23 + EO23I  scal + -EOI          EO3       -EO123          EO1            0          EO2            0            0  E123 + EO123I  -E1 + -EO1I         EO13  -E2 + -EO2I         EO23          -EO        -EO12            0            0            0  -E12 + EO12I       -EO123          EO1          EO2            0         EO12 ;
         EOI          -EI        -EO1I        -EO2I        -EO3I           EO         -E1I         -E2I         -E3I         scal       -EO12I       -EO13I          EO1       -EO23I          EO2          EO3         E12I         E13I          -E1         E23I          -E2          -E3       EO123I        -EO12        -EO13        -EO23        E123I         -E12         -E13         -E23       -EO123         E123 ;
         E12        -E12I          -E2           E1         E123        -EO12          E2I         -E1I       -E123I        EO12I         scal          E23          EO2         -E13         -EO1       -EO123          -EI        -E23I        -EO2I         E13I         EO1I       EO123I           E3          -EO        -EO23         EO13         -E3I          EOI        EO23I       -EO13I         -EO3         EO3I ;
         E13        -E13I          -E3        -E123           E1        -EO13          E3I        E123I         -E1I        EO13I         -E23         scal          EO3          E12        EO123         -EO1         E23I          -EI        -EO3I        -E12I      -EO123I         EO1I          -E2         EO23          -EO        -EO12          E2I       -EO23I          EOI        EO12I          EO2        -EO2I ;
         E1I            0          -EI        -E12I        -E13I   E1 + -EO1I            0            0            0         -E1I         -E2I         -E3I   scal + EOI       -E123I  E12 + EO12I  E13 + EO13I            0            0          -EI            0        -E12I        -E13I         E23I   -E2 + EO2I   -E3 + EO3I  -E123 + EO123I            0          E2I          E3I        E123I  -E23 + -EO23I         E23I ;
         E23        -E23I         E123          -E3           E2        -EO23       -E123I          E3I         -E2I        EO23I          E13         -E12       -EO123         scal          EO3         -EO2        -E13I         E12I       EO123I          -EI        -EO3I         EO2I           E1        -EO13         EO12          -EO         -E1I        EO13I       -EO12I          EOI         -EO1         EO1I ;
         E2I            0         E12I          -EI        -E23I   E2 + -EO2I            0            0            0         -E2I          E1I        E123I  -E12 + -EO12I         -E3I   scal + EOI  E23 + EO23I            0            0         E12I            0          -EI        -E23I        -E13I   E1 + -EO1I  E123 + -EO123I   -E3 + EO3I            0         -E1I       -E123I          E3I  E13 + EO13I        -E13I ;
         E3I            0         E13I         E23I          -EI   E3 + -EO3I            0            0            0         -E3I       -E123I          E1I  -E13 + -EO13I          E2I  -E23 + -EO23I   scal + EOI            0            0         E13I            0         E23I          -EI         E12I  -E123 + EO123I   E1 + -EO1I   E2 + -EO2I            0        E123I         -E1I         -E2I  -E12 + -EO12I         E12I ;
        EO12  E12 + -EO12I         -EO2          EO1        EO123            0    E2 + EO2I  -E1 + -EO1I  -E123 + -EO123I         EO12           EO         EO23            0        -EO13            0            0  scal + -EOI  E23 + -EO23I          EO2  -E13 + EO13I         -EO1       -EO123          EO3            0            0            0  -E3 + -EO3I           EO         EO23        -EO13            0         -EO3 ;
        EO13  E13 + -EO13I         -EO3       -EO123          EO1            0    E3 + EO3I  E123 + EO123I  -E1 + -EO1I         EO13        -EO23           EO            0         EO12            0            0  -E23 + EO23I  scal + -EOI          EO3  E12 + -EO12I        EO123         -EO1         -EO2            0            0            0    E2 + EO2I        -EO23           EO         EO12            0          EO2 ;
        EO1I          E1I         -EOI       -EO12I       -EO13I          EO1           EI         E12I         E13I          -E1        -EO2I        -EO3I           EO      -EO123I         EO12         EO13         -E2I         -E3I         scal       -E123I          E12          E13        EO23I         -EO2         -EO3       -EO123        -E23I           E2           E3         E123        -EO23         -E23 ;
        EO23  E23 + -EO23I        EO123         -EO3          EO2            0  -E123 + -EO123I    E3 + EO3I  -E2 + -EO2I         EO23         EO13        -EO12            0           EO            0            0  E13 + -EO13I  -E12 + EO12I       -EO123  scal + -EOI          EO3         -EO2          EO1            0            0            0  -E1 + -EO1I         EO13        -EO12           EO            0         -EO1 ;
        EO2I          E2I        EO12I         -EOI       -EO23I          EO2        -E12I           EI         E23I          -E2         EO1I       EO123I        -EO12        -EO3I           EO         EO23          E1I        E123I         -E12         -E3I         scal          E23       -EO13I          EO1        EO123         -EO3         E13I          -E1        -E123           E3         EO13          E13 ;
        EO3I          E3I        EO13I        EO23I         -EOI          EO3        -E13I        -E23I           EI          -E3      -EO123I         EO1I        -EO13         EO2I        -EO23           EO       -E123I          E1I         -E13          E2I         -E23         scal        EO12I       -EO123          EO1          EO2        -E12I         E123          -E1          -E2        -EO12         -E12 ;
        E123       -E123I          E23         -E13          E12        EO123        -E23I         E13I        -E12I      -EO123I           E3          -E2         EO23           E1        -EO13         EO12         -E3I          E2I       -EO23I         -E1I        EO13I       -EO12I         scal          EO3         -EO2          EO1          -EI        -EO3I         EO2I        -EO1I           EO         -EOI ;
        E12I            0          E2I         -E1I       -E123I  E12 + EO12I            0            0            0        -E12I           EI         E23I   -E2 + EO2I        -E13I   E1 + -EO1I  E123 + -EO123I            0            0          E2I            0         -E1I       -E123I         -E3I   scal + EOI  E23 + EO23I  -E13 + -EO13I            0          -EI        -E23I         E13I   E3 + -EO3I         -E3I ;
        E13I            0          E3I        E123I         -E1I  E13 + EO13I            0            0            0        -E13I        -E23I           EI   -E3 + EO3I         E12I  -E123 + EO123I   E1 + -EO1I            0            0          E3I            0        E123I         -E1I          E2I  -E23 + -EO23I   scal + EOI  E12 + EO12I            0         E23I          -EI        -E12I   -E2 + EO2I          E2I ;
        E23I            0       -E123I          E3I         -E2I  E23 + EO23I            0            0            0        -E23I         E13I        -E12I  E123 + -EO123I           EI   -E3 + EO3I   E2 + -EO2I            0            0       -E123I            0          E3I         -E2I         -E1I  E13 + EO13I  -E12 + -EO12I   scal + EOI            0        -E13I         E12I          -EI   E1 + -EO1I         -E1I ;
       EO123  -E123 + -EO123I         EO23        -EO13         EO12            0  E23 + -EO23I  -E13 + EO13I  E12 + -EO12I        EO123          EO3         -EO2            0          EO1            0            0  -E3 + -EO3I    E2 + EO2I        -EO23  -E1 + -EO1I         EO13        -EO12           EO            0            0            0  scal + -EOI          EO3         -EO2          EO1            0          -EO ;
       EO12I        -E12I         EO2I        -EO1I      -EO123I         EO12          E2I         -E1I       -E123I          E12          EOI        EO23I         -EO2       -EO13I          EO1        EO123          -EI        -E23I           E2         E13I          -E1        -E123        -EO3I           EO         EO23        -EO13         -E3I         scal          E23         -E13          EO3          -E3 ;
       EO13I        -E13I         EO3I       EO123I        -EO1I         EO13          E3I        E123I         -E1I          E13       -EO23I          EOI         -EO3        EO12I       -EO123          EO1         E23I          -EI           E3        -E12I         E123          -E1         EO2I        -EO23           EO         EO12          E2I         -E23         scal          E12         -EO2           E2 ;
       EO23I        -E23I      -EO123I         EO3I        -EO2I         EO23       -E123I          E3I         -E2I          E23        EO13I       -EO12I        EO123          EOI         -EO3          EO2        -E13I         E12I        -E123          -EI           E3          -E2        -EO1I         EO13        -EO12           EO         -E1I          E13         -E12         scal          EO1          -E1 ;
       E123I            0        -E23I         E13I        -E12I  E123 + -EO123I            0            0            0       -E123I          E3I         -E2I  E23 + EO23I          E1I  -E13 + -EO13I  E12 + EO12I            0            0        -E23I            0         E13I        -E12I          -EI   E3 + -EO3I   -E2 + EO2I   E1 + -EO1I            0         -E3I          E2I         -E1I   scal + EOI          -EI ;
      EO123I        E123I       -EO23I        EO13I       -EO12I        EO123         E23I        -E13I         E12I        -E123         EO3I        -EO2I         EO23         EO1I        -EO13         EO12          E3I         -E2I          E23          E1I         -E13          E12         -EOI          EO3         -EO2          EO1           EI          -E3           E2          -E1           EO         scal ;

      ];
	end
	
        function R = inverse_(A)
            M = productleftexpand_(A);
            if rcond(M) <= eps
                error('Inverse of %s does not exist.', char(A))
            end

           R = CGA(M\CGA(1).m);
        end

        % ***** Norms *****

        function r = norm_(A)
            B = double_(grade_(A.product_(reverse_(A)), 0));
            if B > 0
                r = sqrt(B);
            else
                r = sqrt(-B);
            end
        end

        function r = vnorm_(A)
	    error('The vnorm cannot be performed on CGA elements.');

        end

        function R = normalize_(A)
            if norm_(A) < GA.epsilon_tolerance()
                error("The norm of the element is 0. Cannot normalize.")
            end
            R = A / norm_(A);
        end

        % ***** Equality and Inequality *****

        function b = eq_(A, B)
            b = norm_(A - B) + vnorm_(A - B) < GA.epsilon_tolerance;
            % TODO: Double check that confirming the norm and vnorm are close to 0
            %       is actually sufficient for determining equality.
        end

        function b = eeq_(A, B)
            b = all(A.m == B.m);
        end

        function b = ne_(A, B)
            b = ~eq_(A, B);
        end
        
        % ***** Dual and Reverse*****

        function R = dual_(A)
	    R = -1*leftcontraction(A,I5);
        end

        function R = inversedual_(A)
            error('Inverse dual cannot be performed on CGA elements.');
        end


        function R = reverse_(A)
            R = CGA(  A.m(1), ...
                      [A.m(2); A.m(3); A.m(4); A.m(5); A.m(6)], ...
                    - [A.m(7); A.m(8); A.m(9); A.m(10); A.m(11); A.m(12); A.m(13); A.m(14); A.m(15); A.m(16)], ...
		    - [A.m(17); A.m(18); A.m(19); A.m(20); A.m(21); A.m(22); A.m(23); A.m(24); A.m(25); A.m(26)], ...
                      [A.m(27); A.m(28); A.m(29); A.m(30); A.m(31)], ...
                      A.m(32));
        end

        function R = zeroepsilons_(A)
            R = A;
            for i=1:32
                if abs(R.m(i)) < GA.epsilon_tolerance
                    R.m(i) = 0;
                end
            end
        end

        function R = hodgedual_(A)
            error('Hodge dual cannot be performed on CGA elements.');
        end

        function R = inversehodgedual_(A)
            error('Inverse Hodge dual cannot be performed on CGA elements.');
        end

        function R = jmap_(A)
            error('jmap cannot be performed on CGA elements.');
        end

        function R = poincare_(A)
            error('Poincare cannot be performed on CGA elements.');
        end

        function R = join_(A, B)
            error('join has not been implemented for CGA elements.');
        end

        function R = meet_(A, B)
            R = A.outer_(B);
        end

        function R = conjugate_(A)
            R = CGA(  A.m(1), ...
                    - [A.m(2); A.m(3); A.m(4); A.m(5); A.m(6)], ...
                    - [A.m(7); A.m(8); A.m(9); A.m(10); A.m(11); A.m(12); A.m(13); A.m(14); A.m(15); A.m(16)], ...
		      [A.m(17); A.m(18); A.m(19); A.m(20); A.m(21); A.m(22); A.m(23); A.m(24); A.m(25); A.m(26)], ...
                      [A.m(27); A.m(28); A.m(29); A.m(30); A.m(31)], ...
                    - A.m(32));
        end

        function R = gradeinvolution_(A)
            R = CGA(   A.m(1), ...
                    - [A.m(2); A.m(3); A.m(4); A.m(5); A.m(6)], ...
                      [A.m(7); A.m(8); A.m(9); A.m(10); A.m(11); A.m(12); A.m(13); A.m(14); A.m(15); A.m(16)], ...
		    - [A.m(17); A.m(18); A.m(19); A.m(20); A.m(21); A.m(22); A.m(23); A.m(24); A.m(25); A.m(26)], ...
                      [A.m(27); A.m(28); A.m(29); A.m(30); A.m(31)], ...
                    -  A.m(32));
        end

        function R = grade_(A, n)
            if nargin == 1 || n == -1
                [scalar_nz, vector_nz, bivector_nz, trivector_nz, fourvector_nz, fivevector_nz] = gradestatus_(A);
                if sum([scalar_nz vector_nz bivector_nz trivector_nz fourvector_nz fivevector_nz]) ~= 1
                    R = -1;
                else
                    R = 1*vector_nz + 2*bivector_nz + 3*trivector_nz + 4*fourvector_nz + 5*fivevector_nz; 
                end
            else
                if n == 0
                    R = CGA(A.m(1));
                elseif n == 1
                    R = CGA(0, [A.m(2); A.m(3); A.m(4); A.m(5); A.m(6)], 0, 0, 0, 0);
                elseif n == 2
                    R = CGA(0, 0, [A.m(7); A.m(8); A.m(9); A.m(10); A.m(11); A.m(12); A.m(13); A.m(14); A.m(15); A.m(16)], 0, 0, 0);
                elseif n == 3
                    R = CGA(0, 0, 0, [A.m(17); A.m(18); A.m(19); A.m(20); A.m(21); A.m(22); A.m(23); A.m(24); A.m(25); A.m(26)], 0, 0);
                elseif n == 4
                    R = CGA(0, 0, 0, 0, [A.m(27); A.m(28); A.m(29); A.m(30); A.m(31)], 0);
                elseif n == 5
		            R = CGA(0, 0, 0, 0, 0, A.m(32));
                else
                    R = CGA(0);
                end
            end
        end

        function b = isgrade_(A, g)
            if g == 0
                b = GAisa_(A, "scalar");
            elseif g == 1
                b = GAisa_(A, "vector");
            elseif g == 2
                b = GAisa_(A, "bivector");
            elseif g == 3
                b = GAisa_(A, "trivector");
            elseif g == 4
                b = GAisa_(A, "fourvector");
            elseif g == 5
                b = GAisa_(A, "fivevector");
            elseif g==-1
                b = GAisa_(A, "multivector");
            else
                error('isgrade: invalid grade.');
            end
        end
            
        % TODO: Needs to be upgraded to PGA
        % TODO: make private and wrap for public
        % function r = blade(A)
        %     % blade(A) : return a blade made from the largest portion of a multivector.
        %     A = PGA(A);

        %     s(1) = abs(A.m(1));
        %     s(2) = sqrt(sum(abs(A.m(2:4))));
        %     s(3) = sqrt(sum(abs(A.m(5:7))));
        %     s(4) = abs(A.m(8));
        %     if s(1)>s(2) && s(1)>s(3) && s(1)>s(4)
        %       r = PGA.returnGA_(A.m(1));
        %     elseif s(2)>s(3) && s(2)>s(4)
        %       r = PGA.returnGA_([0; A.m(2); A.m(3); A.m(4); 0; 0; 0; 0]);
        %     elseif s(3)>s(4)
        %       r = PGA.returnGA_([0; 0; 0; 0; A.m(5); A.m(6); A.m(7); 0]);
        %     else
        %       r = PGA.returnGA_([0; 0; 0; 0; 0; 0; 0; A.m(8)]);
        %     end
        % end

        function R = gexp_(A)
            rm = productleftexpand_(A);
            E = expm(rm);
            R = CGA(E(1:32,1));
        end

        function R = glog_(A)
            rm = productleftexpand_(A);
            L = logm(rm);
            R = CGA(L(1:32, 1));
        end

        function R = sqrt_(A)
            rm = productleftexpand_(A);
            S = sqrtm(rm);
            R = CGA(S(1:32, 1));

            % TODO: This is an implementation of equation (90) in PGA4CS.
            %       However, it doesn't resolve the issue with sqrt(PGA(-1)) since we get
            %       0 in the denominator.
            %denom = double_(2*(1 + grade_(A, 0)));
            %R = ((1+A)/sqrt(denom))*(1 - grade_(A, 4)/denom);
        end

        % TODO: Decide behaviour for non-points for get functions.

        function r = getx_(A)
            %GETX_ - A private function for computing the x coordinate of a CGA element.
            %   Returns the x coordinate of a point. Non-points return an error.

            r = A.m(3)/A.m(2);
            % 3 is the position of e1
            % 2 is the position of no
        end
        
        function r = gety_(A)
            %GETY_ - A private function for computing the y coordinate of a CGA element.
            %   Returns the y coordinate of a point. Non-points return an error.

            r = A.m(4)/A.m(2);
            % 4 is the position of e2
            % 2 is the position of no
        end
        
        function r = getz_(A)
            %GETZ_ - A private function for computing the z coordinate of a CGA element.
            %   Returns the z coordinate of a point. Non-points return an error.

            r = A.m(5)/A.m(2);
            % 5 is the position of e3
            % 2 is the position of no
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

            [s, pl] = GA.charifyval_(p.m(2), 'no', s, pl);
            [s, pl] = GA.charifyval_(p.m(3), 'e1', s, pl);
            [s, pl] = GA.charifyval_(p.m(4), 'e2', s, pl);
            [s, pl] = GA.charifyval_(p.m(5), 'e3', s, pl);
            [s, pl] = GA.charifyval_(p.m(6), 'ni', s, pl);

	        [s, pl] = GA.charifyval_(p.m(7), 'no^e1', s, pl);
            [s, pl] = GA.charifyval_(p.m(8), 'no^e2', s, pl);
	        [s, pl] = GA.charifyval_(p.m(9), 'no^e3', s, pl);
 	        [s, pl] = GA.charifyval_(p.m(10), 'no^ni', s, pl);
            if GA.compact_notation()
                if ~PGA.increasing_order()
                    [s, pl] = GA.charifyval_(p.m(11), 'e23', s, pl);
                    [s, pl] = GA.charifyval_(-p.m(10), 'e31', s, pl);
                    [s, pl] = GA.charifyval_(p.m(9), 'e12', s, pl);

                    [s, pl] = GA.charifyval_(-p.m(14), 'e032', s, pl);
                    [s, pl] = GA.charifyval_(p.m(13), 'e013', s, pl);
                    [s, pl] = GA.charifyval_(-p.m(12), 'e021', s, pl);
                    [s, pl] = GA.charifyval_(p.m(15), 'e123', s, pl);
                else 
                    [s, pl] = GA.charifyval_(p.m(11), 'e12', s, pl);
                    [s, pl] = GA.charifyval_(p.m(12), 'e13', s, pl);
                    [s, pl] = GA.charifyval_(p.m(13), 'e1^ni', s, pl);

                    [s, pl] = GA.charifyval_(p.m(12), 'e012', s, pl);
                    [s, pl] = GA.charifyval_(p.m(13), 'e013', s, pl);
                    [s, pl] = GA.charifyval_(p.m(14), 'e023', s, pl);
                    [s, pl] = GA.charifyval_(p.m(15), 'e123', s, pl);
                end

                if GA.compact_pseudoscalar()
                    [s, pl] = GA.charifyval_(p.m(16), 'I5', s, pl);
                else
                    [s, pl] = GA.charifyval_(p.m(16), 'no^e123^ni', s, pl);
                end
            else
                if ~PGA.increasing_order()
                    [s, pl] = GA.charifyval_(p.m(14), 'e2^e3', s, pl);
                    [s, pl] = GA.charifyval_(-p.m(12), 'e3^e1', s, pl);
                    [s, pl] = GA.charifyval_(p.m(11), 'e1^e2', s, pl);
                    [s, pl] = GA.charifyval_(p.m(13), 'e1^ni', s, pl);
                    [s, pl] = GA.charifyval_(p.m(15), 'e2^ni', s, pl);
                    [s, pl] = GA.charifyval_(p.m(16), 'e3^ni', s, pl);
		    
                    [s, pl] = GA.charifyval_(p.m(20), 'no^e2^e3', s, pl);
                    [s, pl] = GA.charifyval_(-p.m(18), 'no^e3^e1', s, pl);
                    [s, pl] = GA.charifyval_(p.m(17), 'no^e1^e2', s, pl);
                    [s, pl] = GA.charifyval_(p.m(19), 'no^e1^ni', s, pl);
                    [s, pl] = GA.charifyval_(p.m(21), 'no^e2^ni', s, pl);
                    [s, pl] = GA.charifyval_(p.m(22), 'no^e3^ni', s, pl);
                    [s, pl] = GA.charifyval_(p.m(23), 'e1^e2^e3', s, pl);
                    [s, pl] = GA.charifyval_(p.m(26), 'e2^e3^ni', s, pl);
                    [s, pl] = GA.charifyval_(-p.m(25), 'e3^e1^ni', s, pl);
                    [s, pl] = GA.charifyval_(p.m(24), 'e1^e2^ni', s, pl);

                    [s, pl] = GA.charifyval_(p.m(30), 'no^e2^e3^ni', s, pl);
                    [s, pl] = GA.charifyval_(-p.m(29), 'no^e3^e1^ni', s, pl);
                    [s, pl] = GA.charifyval_(p.m(28), 'no^e1^e2^ni', s, pl);
                    [s, pl] = GA.charifyval_(p.m(27), 'no^e1^e2^e3', s, pl);
                    [s, pl] = GA.charifyval_(p.m(31), 'e1^e2^e3^ni', s, pl);
                else
                    [s, pl] = GA.charifyval_(p.m(11), 'e1^e2', s, pl);
                    [s, pl] = GA.charifyval_(p.m(12), 'e1^e3', s, pl);
                    [s, pl] = GA.charifyval_(p.m(13), 'e1^ni', s, pl);
                    [s, pl] = GA.charifyval_(p.m(14), 'e2^e3', s, pl);
                    [s, pl] = GA.charifyval_(p.m(15), 'e2^ni', s, pl);
                    [s, pl] = GA.charifyval_(p.m(16), 'e3^ni', s, pl);

                    [s, pl] = GA.charifyval_(p.m(17), 'no^e1^e2', s, pl);
                    [s, pl] = GA.charifyval_(p.m(18), 'no^e1^e3', s, pl);
                    [s, pl] = GA.charifyval_(p.m(19), 'no^e1^ni', s, pl);
                    [s, pl] = GA.charifyval_(p.m(20), 'no^e2^e3', s, pl);
                    [s, pl] = GA.charifyval_(p.m(21), 'no^e2^ni', s, pl);
                    [s, pl] = GA.charifyval_(p.m(22), 'no^e3^ni', s, pl);
                    [s, pl] = GA.charifyval_(p.m(23), 'e1^e2^e3', s, pl);
                    [s, pl] = GA.charifyval_(p.m(24), 'e1^e2^ni', s, pl);
                    [s, pl] = GA.charifyval_(p.m(25), 'e1^e3^ni', s, pl);
                    [s, pl] = GA.charifyval_(p.m(26), 'e2^e3^ni', s, pl);

                    [s, pl] = GA.charifyval_(p.m(27), 'no^e1^e2^e3', s, pl);
                    [s, pl] = GA.charifyval_(p.m(28), 'no^e1^e2^ni', s, pl);
                    [s, pl] = GA.charifyval_(p.m(29), 'no^e1^e3^ni', s, pl);
                    [s, pl] = GA.charifyval_(p.m(30), 'no^e2^e3^ni', s, pl);
                    [s, pl] = GA.charifyval_(p.m(31), 'e1^e2^e3^ni', s, pl);
                end

                if GA.compact_pseudoscalar()
                    [s, pl] = GA.charifyval_(p.m(32), 'I5', s, pl);
                else
                    [s, pl] = GA.charifyval_(p.m(32), 'no^e1^e2^e3^ni', s, pl);
                end
            end
            
            if strcmp(pl, ' ')
                s = '     0';
            end
        end
    end

    % ******************** Public Static Methods ********************
    methods (Access = public, Static)
        function e = elements()
            if PGA.increasing_order()
                e = [PGA(1), e0(PGA), e1(PGA), e2(PGA), e3(PGA), e01(PGA), e02(PGA), e03(PGA),...
                     e12(PGA), e13(PGA), e23(PGA), e012(PGA), e013(PGA), e023(PGA), e123(PGA), e0123(PGA)];
            else 
                e = [PGA(1), e0(PGA), e1(PGA), e2(PGA), e3(PGA), e01(PGA), e02(PGA), e03(PGA),...
                     e23(PGA), e31(PGA), e12(PGA), e032(PGA), e013(PGA), e021(PGA), e123(PGA), e0123(PGA)];
            end
        end

        function s = modelname()
            s = "CGA";
        end

        function R = cast(A)
            if isa(A, 'CGA')
                R = A;
            elseif isa(A, 'double')
                if GA.autoscalar
                    R = CGA(A);
                else 
                    error('Implicit conversion between a double and CGA is disabled. Run "help GA.autoscalar" for more information.')
                end
            else
                error(['Cannot implictly convert from ' class(A) ' to CGA'])
            end
        end

        function r = getzero()
            r = CGA(0);
        end
    end


    % ******************** Public Methods ********************

    methods (Access = public)
        function obj = CGA(m0, m1, m2, m3, m4, m5)
            %CGA - The constructor for PGA elements.
            %   If no arugment is provided, the 0 element in PGA is returned.
            %   If 1 arugment is provided and it is a PGA element, the element itself will
            %   be returned. If the argument is a double, a PGA element of that scalar will
            %   be returned.
            %   If 6 arguments are provided, it is assumed they have the following format:
            %   The first argument is a double for the scalar portion of the multivector
            %   The second argument is [c0, c1, c2, c3, c4], where ci is a double corresponding
            %   to the coefficient for ei.
            %   The third argument is [c01, c02, c03, c12, c13, c23], where cij is a double
            %   corresponding to the coefficient for eij.
            %   The fourth argument is [c012, c013, c023, c123], where cijk is a double
            %   corresponding to the coefficient for eijk.
            %   The fifth argument is a double corresponding to the coefficient for e0123.
            %   For any of the arguments, one can optionally simply put 0 to have zero
            %   for all the corresponding coefficients. 
            %
            %   It is not generally recommended to construct PGA elements this way.
            %   Instead, consider writing equations of the form
            %                              e1 + e2^e3
            %   while in the PGA model (see help GA.settings and help GA.model)
            %   or while not in the PGA model as
            %                              e1(PGA) + e2(PGA)^e3(PGA)
            if nargin == 0
                obj = CGA(0);
            elseif nargin == 6
                if m1 == 0
                    m1 = zeros(5, 1);
                end
                if m2 == 0
                    m2 = zeros(10, 1);
                end
                if m3 == 0
                    m3 = zeros(10, 1);
		end
		if m4 == 0
                    m4 = zeros(5, 1);
                end
                obj.m = [m0; 
                         m1(1); m1(2); m1(3); m1(4); m1(5);
                         m2(1); m2(2); m2(3); m2(4); m2(5); m2(6); m2(7); m2(8); m2(9); m2(10); 
                         m3(1); m3(2); m3(3); m3(4); m3(5); m3(6); m3(7); m3(8); m3(9); m3(10); 
                         m4(1); m4(2); m4(3); m4(4); m4(5);
                         m5];
            elseif nargin == 1
                if isa(m0, 'CGA')
                    obj = m0;
                elseif isa(m0, 'double')
                    if size(m0, 1) == 1 & size(m0, 2) == 1
                        % User has provided a scalar
                        obj.m = [m0; 
                                    0; 0; 0; 0; 0;
                                    0; 0; 0; 0; 0; 0; 0; 0; 0; 0;
                                    0; 0; 0; 0; 0; 0; 0; 0; 0; 0;
                                    0; 0; 0; 0; 0;
                                    0];
                    elseif size(m0, 1) == 1 & size(m0, 2) == 32
                        % User has provided a column vector
                        obj.m = m0';
                    elseif size(m0, 1) == 32 & size(m0, 2) == 1
                        % User has prodived a row vector
                        obj.m = m0;
                    else
                        error('Bad CGA argument: Unexpected array size.\nExpected size is either 1x1, 32x1 or 1x32.\nCurrent size is: %dx%d', size(m0, 1), size(m0, 2))
                    end
                else
                    error('Bad CGA argument: Invalid input type. Class is currently: %s.', class(m0))
                end
            else 
                error('Bad CGA argument: Invalid number of arguments.')
            end
        end

        function rm = matrix(A)
            rm = A.m;
        end

        function b = GAisa(A, t)
            %GAISA - Determines in a multivector and a string representing a type of multivector
            %   and returns true if the multivector is of that type.
            %   In CGA, valid types are:
            %   scalar, vector, plane, bivector, line, trivector, point, quadvector, quintvector,
            %   pseudoscalar, multivector

            arguments
                A CGA;
                t (1, 1) string;
            end
            
            b = GAisa_(A, t);
        end

        function R = CGAcast(A)
            R = A;
        end

        function R = PGAcast(A)
	        display('PGAcast: Still to be written\n');
            R = A;
        end

        function R = OGAcast(A)
            [scal,  EO, E1, E2, E3, EI,  EO1, EO2, EO3, EOI, E12, E13, E1I, E23, E2I, E3I,  EO12, EO13, EO1I, EO23, EO2I, EO3I, E123, E12I, E13I, E23I, EO123, EO12I, EO13I, EO23I, E123I, EO123I] = decompose_(A);
            R = OGA(scal, [E1, E2, E3], [E12, E13, E23], E123);
        end

        function drawTP(A, center, c)
            arguments
                A PGA;
		        center PGA;
                c = [];
            end
            if GAisa(A, 'plane')
                if isempty(c)
                    % TODO: make default colour of plane changable.
                    c = 'g';
                end

		        dist = double(outer(A,center).noneuclidean.*I3);
                % Check if center on plane A; if not, then project onto A
                if abs(dist) > eps
                    nrm = norm(A.euclidean)
                    t = A.euclidean*dist/(nrm*nrm);
                    trans = 1-e0*(t)/2;
                    transi = 1+e0*(t)/2;
                    center = trans*center*transi;
                end
                h = PGABLEDraw.pointingplaneC(A, center, c);
                GAScene.addstillitem(GASceneStillItem(A, h));
            else
                error('Error is not a point, line, or plane. PGA cannot draw it.');
            end
	    end

        function draw(A, varargin)
            arguments
                A CGA;
            end
            arguments (Repeating)
                varargin
            end
            GAScene.usefigure();
	    
            A = zeroepsilons_(A);
            
            [scal,  EO, E1, E2, E3, EI,  EO1, EO2, EO3, EOI, E12, E13, E1I, E23, E2I, E3I,  EO12, EO13, EO1I, EO23, EO2I, EO3I, E123, E12I, E13I, E23I, EO123, EO12I, EO13I, EO23I, E123I, EO123I] = decompose_(A);
            if GAisa(A, 'point')

                % Custom input handling
                argsize = size(varargin, 2);
                if argsize == 1
                    if isa(varargin{1}, "char")
                        varargin = ['FaceColor', varargin];
                    end
                end

                updated_varargin = PGABLEDraw.defaultvarargin('FaceColor', 'y', varargin{:});
%                if euclidean(A) == 0
                if 0
                    % TODO perhaps make drawing as part of the constructor of the dynamic item
                    h = PGA.drawvanishingpoint(A, updated_varargin{:});
                    GAScene.adddynamicitem(GASceneDynamicItem(A, h, @()PGA.drawvanishingpoint(A, updated_varargin{:})));
                else
                    h = PGABLEDraw.octahedron(A, CGA.pointsize(), updated_varargin{:});
                    GAScene.addstillitem(GASceneStillItem(A, h));
                end
                
            elseif GAisa(A, 'sphere')
                % don't worry about dual vs primal spheres for now
                if isgrade(A,4)
                    A = dual(A);
                end
                % Only worry about sign of weight (no) for now
                if A.nocoeff() > 0
                    outn = 1;
                else
                    outn = 0;
                end
                A = (1./A.nocoeff())*A;
                cx = A.getx(); cy = A.gety(); cz = A.getz();
                len = sqrt(cx*cx + cy*cy + cz*cz);
                r = 2*A.nicoeff()-len;
                r = sqrt(double(1./double((A^ni)*(A^ni)) * (A*A)));
                % Don't worry about imaginary vs real spheres for now
                if r<0
                    r = -1*r;
                end
                %r = sqrt(r)
                h = PGABLEDraw.wfsphere(A, r, outn, varargin{:});
                GAScene.addstillitem(GASceneStillItem(A,h));
            elseif GAisa(A, 'circle')
                if grade(A) == 2
                    A = dual(A);
                    [scal,  EO, E1, E2, E3, EI,  EO1, EO2, EO3, EOI, E12, E13, E1I, E23, E2I, E3I,  EO12, EO13, EO1I, EO23, EO2I, EO3I, E123, E12I, E13I, E23I, EO123, EO12I, EO13I, EO23I, E123I, EO123I] = decompose_(A);
                end
                %A = (1./norm(A))*A


	    argsize = size(varargin, 2);
	    if argsize == 1
	       if isa(varargin{1}, "char")
	         varargin = ['Color', varargin];
	       end
	    end
	    updated_varargin = PGABLEDraw.defaultvarargin('Color', 'y', varargin{:});
	    
                nx = EO23; ny = -EO13; nz = EO12;

                nlen = sqrt(nx*nx + ny*ny + nz*nz);
                unx = nx/nlen; uny = ny/nlen; unz = nz/nlen;
                nA = 1./nlen*A;

                % for testing
                %pln = A^ni
                %draw(2*pln)

                cp = A*ni*A;
                cp = (1./cp.nocoeff())*cp;
                cpx = cp.e1coeff(); cpy = cp.e2coeff(); cpz = cp.e3coeff();
                cpsq = cpx*cpx + cpy*cpy + cpz*cpz;

 		%r = sqrt(2*E1I/nx-cpsq);

                %r = sqrt(4*(cpsq-(E23I/nx-E13I/ny+E12I/nz)))
                %r = sqrt(cpsq-2*E23I/nx)
                %r = sqrt(cpsq-2*E12I/nz)
        %		%r = sqrt(-1/sqrt(norm(A))*double((dual(A)*dual(A))))
                % Works for a sphere...
                r = sqrt(-1*double(1./double((A^ni)*(A^ni)) * (A*A)));
		isImaginary=0;
		if ( imag(r) ~= 0 )
		  r = imag(r);
		  isImaginary=1
		end

                if abs(unx)<=abs(uny) && abs(unx)<=abs(unz)
                    vvec = cross([1 0 0],[unx,uny,unz]);
                    wvec = cross(vvec,[unx,uny,unz]);
                elseif abs(uny) <= abs(unz) && abs(uny) <= abs(unx)
                    vvec = cross([0 1 0],[unx,uny,unz]);
                    wvec = cross(vvec,[unx,uny,unz]);
                else
                    vvec = cross([0 0 1],[unx,uny,unz]);
                    wvec = cross(vvec,[unx,uny,unz]);
                end
                vvec = 1./norm(vvec)*vvec;
                wvec = 1./norm(wvec)*wvec;
                %draw(cp) % for testing
                ptB = r*vvec;
                for ii=0:51
                    ang = ii/50*2*pi;
                    ptA = r*cos(ang)*vvec + r*sin(ang)*wvec;
                    if ii ~= 0 && (~isImaginary || mod(ii,2)==0)
                        plot3([ptB(1)+cpx ptA(1)+cpx], ...
                            [ptB(2)+cpy ptA(2)+cpy], ...
                            [ptB(3)+cpz ptA(3)+cpz], 'LineWidth', 2, updated_varargin{:});
                        hold on;
                    end
                    ptB = ptA;
                end
            elseif GAisa(A, 'line')

                offset = [];
                % Custom input handling
                argsize = size(varargin, 2);
                if argsize == 1
                    if isa(varargin{1}, "char")
                        varargin = ['Color', varargin];
                    elseif isa(varargin{1}, "GA")
                        offset = varargin{1};
                        varargin = {};
                    end
                elseif argsize == 2
                    if isa(varargin{1}, "GA") && isa(varargin{2}, "char")
                        offset = varargin{1};
                        varargin{1} = 'Color';
                    elseif isa(varargin{1}, "char") && isa(varargin{2}, "GA")
                        error("Arguments are in an incorrect order. It should be draw(ELEMENT, OFFSET, COLOR).")
                    end
                end

                updated_varargin = PGABLEDraw.defaultvarargin('Color', 'b', varargin{:});
                updated_varargin = PGABLEDraw.defaultvarargin('LineWidth', 1.5, updated_varargin{:});
    %                if euclidean(A) == 0
    %                    if ~isempty(offset)
    %                        error("Cannot offset lines at infinity. Do not provide an offset argument.")
    %                    end
    %                    % TODO perhaps make drawing as part of the constructor of the dynamic item
    %                    %TODO: make dashedness work
    %                    updated_varargin = [{'--'}, updated_varargin];
    %                    h = PGA.drawvanishingline(A, updated_varargin{:});
    %                    GAScene.adddynamicitem(GASceneDynamicItem(A, h, @()PGA.drawvanishingline(A, updated_varargin{:})));
    %                else

                % Extract direction, vector to origin
                vx = EO1I; vy = EO2I; vz = EO3I;
                mx = E23I; my = -E13I; mz = E12I; % not quite sure of sign
                %len = sqrt(vx*vx+vy*vy+vz*vz)
                %mx = mx/len; my=my/len; mz=mz/len;
                %mx,my,mz
                pgaPc = cross([vx,vy,vz],[mx,my,mz]);
                pgaPc = pgaPc/(vx*vx+vy*vy+vz*vz);
                pgaP = gapoint(pgaPc(1),pgaPc(2),pgaPc(3),PGA);
                pgaV = vx*e1(PGA) + vy*e2(PGA) + vz*e3(PGA);
                pgaA = pgaP .*  pgaV;
                if isempty(offset)
                    h = PGABLEDraw.hairyline(pgaA, updated_varargin{:});
                else 
                    h = PGABLEDraw.hairylineC(pgaA, offset, updated_varargin{:});
                end
                        
                GAScene.addstillitem(GASceneStillItem(pgaA, h));
    %                end

            elseif GAisa(A, 'plane')

                % Don't worry about dual planes for now
                if isgrade(A,4)
                    A = dual(A);
                end
                    % Calculating center
    %                plane = normalize(A);
    %                if euclidean(plane) == 0
    %                    error("This is a plane at infinity. Cannot currently display this object.");
    %                end
    %                delta = -e0coeff(A);
                support_vec = euclidean(A) ;
                delta = A.nicoeff();
                center = delta*support_vec + no(CGA);

                % Setting the default offset to be the desired center
                offset = center;

                % Custom input handling
                argsize = size(varargin, 2);
                if argsize == 1
                    if isa(varargin{1}, "char")
                        varargin = ['FaceColor', varargin];
                    elseif isa(varargin{1}, "GA")
                        offset = varargin{1};
                        varargin = {};
                    end
                elseif argsize == 2
                    if isa(varargin{1}, "GA") && isa(varargin{2}, "char")
                        offset = varargin{1};
                        varargin{1} = 'FaceColor';
                    elseif isa(varargin{1}, "char") && isa(varargin{2}, "GA")
                        error("Arguments are in an incorrect order. It should be draw(ELEMENT, OFFSET, COLOR).")
                    end
                end


                updated_varargin = PGABLEDraw.defaultvarargin('FaceColor', 'g', varargin{:});
                updated_varargin = PGABLEDraw.defaultvarargin('FaceAlpha', 0.5, updated_varargin{:});
                % Convert CGA to PGA
                pgaoffset = gapoint(offset.getx(), offset.gety(), offset.getz(), PGA);
                pgaA = support_vec.e1coeff()*e1(PGA) +...
                support_vec.e2coeff()*e2(PGA) +...
                support_vec.e3coeff()*e3(PGA) -...
                delta*e0(PGA)
                h = PGABLEDraw.pointingplaneC(pgaA, pgaoffset, updated_varargin{:});
                GAScene.addstillitem(GASceneStillItem(A, h));
            else
                    error('This element is not of a type that can be drawn in CGAs.');
            end
        end

        function R = euclidean(A)
            %EUCLIDEAN - Returns the euclidean portion of the multivector.
	        % For CGA, this is the e1, e2, e3 portion

            [scal,  EO, E1, E2, E3, EI,  EO1, EO2, EO3, EOI, E12, E13, E1I, E23, E2I, E3I,  EO12, EO13, EO1I, EO23, EO2I, EO3I, E123, E12I, E13I, E23I, EO123, EO12I, EO13I, EO23I, E123I, EO123I] = decompose_(A);
            R = CGA(0, [0, E1, E2, E3,0 ], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0], 0); 
        end

        function R = noneuclidean(A)
            %NONEUCLIDEAN - Returns the non-euclidean portion of the multivector.
	        % For CGA, this is the no, ni portion

            [scal,  EO, E1, E2, E3, EI,  EO1, EO2, EO3, EOI, E12, E13, E1I, E23, E2I, E3I,  EO12, EO13, EO1I, EO23, EO2I, EO3I, E123, E12I, E13I, E23I, EO123, EO12I, EO13I, EO23I, E123I, EO123I] = decompose_(A);
            R = CGA(0, [EO, 0, 0, 0, EI], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 0, 0, 0);
        end

        function r = nocoeff(A)
            %NOCOEFF - Returns the coefficient of no.

            M = matrix(A);
            r = M(2);
        end

        function r = e1coeff(A)
            %E1COEFF - Returns the coefficient of e1.

            M = matrix(A);
            r = M(3);
        end

        function r = e2coeff(A)
            %E2COEFF - Returns the coefficient of e2.

            M = matrix(A);
            r = M(4);
        end

        function r = e3coeff(A)
            %E3COEFF - Returns the coefficient of e3.

            M = matrix(A);
            r = M(5);
        end

        function r = nicoeff(A)
            %NICOEFF - Returns the coefficient of ni.

            M = matrix(A);
            r = M(6);
        end

        function r = noe123nicoeff(A)
            %NOE123NICOEFF - Returns the coefficient of noe123ni.

            M = matrix(A);
            r = M(32);
        end
    end
end

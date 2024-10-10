classdef GAScene
    %GASCENE - A class which renders of a figure for a specific set of GA objects.
    %   Only the following methods should be used by a beginner:
    %      • GAScene.clearitems()      to clear items in the GA Scene
    %      • GAScene.displayitems()    to see a list of items displayed in the GAScene
    %   If you need to make modifications to the figure itself, you can acquire the figure handle by running
    %      • GAScene.getfigure()
    %
    % See also GA, OGA, PGA.

    % PGABLE, Copyright (c) 2024, University of Waterloo
    % Copying, use and development for non-commercial purposes permitted.
    %          All rights for commercial use reserved; for more information
    %          contact Stephen Mann (smann@uwaterloo.ca)
    %
    %          This software is unsupported.

    properties (Access = private)
        fig;   % The figure to be drawn to
        items; % The items that exist in the scene
        dynamic_items; % The items that need updating after moving the scene
        scenename;
    end

    % ******************** Public Static Methods ********************

    methods (Access = public, Static)

        
        %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%
        %     General Draw Methods     %
        %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%

        function h = drawarrow(p1, p2, c)

            % Drawing an arrow from point p1 to point p2 in colour c
            arguments
                p1 PGA;
                p2 PGA;
                c = 'b';
            end

            % TODO: Eventually use PGA to perform rendering.
            % Sadly, for now we have reverting to rendering by hand.

            GAScene.usefigure();
            hold on

            p1x = p1.getx();
            p1y = p1.gety();
            p1z = p1.getz();

            p2x = p2.getx();
            p2y = p2.gety();
            p2z = p2.getz();
            
            p = (p2x - p1x)*e1(OGA) + (p2y - p1y)*e2(OGA) + (p2z - p1z)*e3(OGA);
            biv = e1(OGA)^p;
            if biv == 0
                dir1 = normalize(dual(e2(OGA)^p));
            else
                dir1 = normalize(dual(biv));
            end

            beta = 1/16;
            % Setting the size of dir1 to be proportional to p
            dir1 = beta*norm(p)*dir1;

            dir1x = dir1.getx();
            dir1y = dir1.gety();
            dir1z = dir1.getz();

            dir2 = normalize(dual(dir1^p));
            dir2 = beta*norm(p)*dir2;

            dir2x = dir2.getx();
            dir2y = dir2.gety();
            dir2z = dir2.getz();


            % DRAWING
            % First, plot line from point 1 to point 2
            h = GAScene.plotline({p1 p2}, c);

            % Then, for each angle, we want to rotate the point around dir1 and dir2
            phi = 0.75; % The base of the arrow head is 0.75 of the way
            abx = (1-phi) * p1x + phi * p2x;
            aby = (1-phi) * p1y + phi * p2y;
            abz = (1-phi) * p1z + phi * p2z;

            g = 8;
            theta = 2*pi/g;
            for i=1:g
                angle = theta*i;
                ht = plot3([p2x, abx + sin(angle)*dir1x + cos(angle)*dir2x], ...
                           [p2y, aby + sin(angle)*dir1y + cos(angle)*dir2y], ...
                           [p2z, abz + sin(angle)*dir1z + cos(angle)*dir2z], ...
                        c);
                h = [h ht];
            end

            for i=0:g
                angle1 = theta*i;
                angle2 = angle1+theta;
                ht = plot3([abx + sin(angle1)*dir1x + cos(angle1)*dir2x, abx + sin(angle2)*dir1x + cos(angle2)*dir2x], ...
                           [aby + sin(angle1)*dir1y + cos(angle1)*dir2y, aby + sin(angle2)*dir1y + cos(angle2)*dir2y], ...
                           [abz + sin(angle1)*dir1z + cos(angle1)*dir2z, abz + sin(angle2)*dir1z + cos(angle2)*dir2z], ...
                        c);
                h = [h ht];
            end

        end

       
        function polygon_(Poly, BV, position, c)

            % The polygon struct assumes the polygon is in the (x,y)-plane. 
            % Polygon is a struct containing the information of the polygon. It contains:
            %   - X, Y, Z vectors which are list of (x, y, z) coordinates for each point in the polygon
            %   - CX and CY, which is the (x, y) position of the center of the polygon

            % TODO: Consider making a public version of this method.
        
            translate = position/origin(PGA);

            rotate = e12(PGA) / BV;
            rotate = rotate/norm(rotate);
            rotate = PGA(sqrt(PGA(rotate)));
            
            for i=1:length(Poly.X)
                p = hodgedual(((Poly.X(i)-Poly.CX)*e1(PGA)*0.01 + (Poly.Y(i)-Poly.CY)*e2(PGA)*0.01) + e0(PGA));
                p = rotate * p * inverse(rotate);
                p = translate * p;
                
                pts{i} = p;
            end
            
            GAScene.patch(pts, c);
        end

        % TODO: Consider making this a public method. Perhaps rename to not be confused with MATLAB's patch function.
        %       Or, potentially pass alpha and line_style in the same style as Matlab to ensure compatability.
        % TODO: perhaps add underscore if not making public
        function handle = patch(pts, c, alpha, line_style)
            arguments
                pts;
                c = 'y';
                alpha = 1;
                line_style = '-';
            end

            for i = 1:length(pts)
                p = pts{i};
                % TODO: preallocate space of x, y, z for speed.
                x(i) = p.getx();
                y(i) = p.gety();
                z(i) = p.getz();
            end
            handle = patch(x, y, z, c, 'FaceAlpha', alpha, 'LineStyle', line_style);
        end

        function h = plotline(A, c, isdashed, isthick)
            % plot3(A, c): plot a 3d line connecting end points of A in color c
            arguments
                A;
                c = [0, 0, 1, 1];
                isdashed = false;
                isthick = false;
            end

            if length(A) == 0;
                h = [];
                return;
            end
            
            for i = 1:length(A)
                p = A{i};
                x(i) = p.getx();
                y(i) = p.gety();
                z(i) = p.getz();
            end

            if isdashed
                if isthick
                    h = plot3(x, y, z, '--', 'LineWidth', 1.5);
                else
                    h = plot3(x, y, z, '--');
                end
            else
                if isthick
                    h = plot3(x, y, z, 'LineWidth', 1.5);
                else
                    h = plot3(x, y, z);
                end
            end
            h.Color = c;
        end

        %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%
        %     Generic Draw Methods     %
        %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%
        function h = drawoctahedron(center_point, radius, c)
            arguments 
                center_point PGA;
                radius double;
                c;
            end

            GAScene.usefigure();
            hold on

            p_e1m = (point(-radius, 0, 0, PGA)/origin(PGA))*center_point;
            p_e1p = (point(radius, 0, 0, PGA)/origin(PGA))*center_point;
            p_e2m = (point(0, -radius, 0, PGA)/origin(PGA))*center_point;
            p_e2p = (point(0, radius, 0, PGA)/origin(PGA))*center_point;
            p_e3m = (point(0, 0, -radius, PGA)/origin(PGA))*center_point;
            p_e3p = (point(0, 0, radius, PGA)/origin(PGA))*center_point;

            h = [GAScene.patch({p_e1m, p_e2m, p_e3m}, c), ...
                 GAScene.patch({p_e1m, p_e2m, p_e3p}, c), ...
                 GAScene.patch({p_e1m, p_e2p, p_e3m}, c), ...
                 GAScene.patch({p_e1m, p_e2p, p_e3p}, c), ...
                 GAScene.patch({p_e1p, p_e2m, p_e3m}, c), ...
                 GAScene.patch({p_e1p, p_e2m, p_e3p}, c), ...
                 GAScene.patch({p_e1p, p_e2p, p_e3m}, c), ...
                 GAScene.patch({p_e1p, p_e2p, p_e3p}, c)];
        end

        function h = drawhairyline(line, c)

            GAScene.usefigure();
            hold on

            % TODO: Currently normalizing the line. We should perhaps shouldn't need to do this.
            line = normalize(line);
            % Direction the line is pointing
            dir = euclidean(line)/I3(PGA);
            % Closest point to the origin
            p = line/dir;
            % Translation in direction from origin
            translation = ihd(dir + e0(PGA))/origin(PGA);
            trans = sqrt(translation);

            point_a = trans * p * inverse(trans);
            point_b = inverse(trans) * p * trans;
            
            
            h = GAScene.plotline({point_a, point_b}, c, false, true);

            % Drawing hairs
            num_hairs = 20;

            hair_trans = gexp(2*glog(trans)/num_hairs);

            phi = pi/3;
            % TODO: PGA4CS says to use -phi. It is not clear to me why this should be the case.
            hair_rot = gexp(-phi*normalize(line)/2);

            hair_trans_rot = hair_trans * hair_rot;
            

            % Creating hair base and tip.
            hair_base = point_b;
            hair_tip = point(0.05, 0, 0, PGA);
            
            % Rotate
            hair_tip_rot = sqrt(euclidean(line)/(e1(PGA)^e2(PGA)));
            hair_tip = hair_tip_rot * hair_tip * inverse(hair_tip_rot);

            % Translate
            hair_tip_trans = sqrt(point_b/origin(PGA));
            hair_tip = hair_tip_trans * hair_tip * inverse(hair_tip_trans);

            for i = 1:num_hairs
                hair_base = hair_trans_rot * hair_base * inverse(hair_trans_rot);
                h = [h GAScene.plotline({hair_base, hair_tip}, c, false, true)];
                hair_tip = hair_trans_rot * hair_tip * inverse(hair_trans_rot);
            end
        end

        function h = drawpointingplane(plane, c)
            % TODO: verify input is PGA vector
            arguments
                plane;
                c = 'g';
            end

            GAScene.usefigure();
            hold on

            % TODO: Can currently only draw normalized planes. Perhaps should visualize non-unitness somehow.
            beta = sqrt(norm(plane));
            plane = normalize(plane);

            if euclidean(plane) == 0
                error("This is a plane at infinity. Cannot currently display this object.");
            end


            delta = -e0coeff(plane);

            support_vec = euclidean(plane); 

            center = ihd(delta*support_vec + e0(PGA));

            sv = OGAcast(support_vec);
            
            biv = e1(OGA)^sv;
            if biv == 0
                dir1 = normalize(dual(e2(OGA)^sv));
            else
                dir1 = normalize(dual(biv));
            end

            dir2 = normalize(dual(dir1^sv));

            dir1 = beta*dir1;
            dir2 = beta*dir2;

            cx = center.getx();
            cy = center.gety();
            cz = center.getz();

            dir1x = dir1.getx();
            dir1y = dir1.gety();
            dir1z = dir1.getz();

            dir2x = dir2.getx();
            dir2y = dir2.gety();
            dir2z = dir2.getz();

            p1t = point(cx + dir1x + dir2x, cy + dir1y + dir2y, cz + dir1z + dir2z, PGA);
            p2t = point(cx - dir1x + dir2x, cy - dir1y + dir2y, cz - dir1z + dir2z, PGA);
            p3t = point(cx - dir1x - dir2x, cy - dir1y - dir2y, cz - dir1z - dir2z, PGA);
            p4t = point(cx + dir1x - dir2x, cy + dir1y - dir2y, cz + dir1z - dir2z, PGA);
            
            h = GAScene.patch({p1t, p2t, p3t, p4t}, c, 0.5);

            
            

            % gs is the grid size. 1 -> one arrow on each corner
            %                      2 -> arrows in 3x3 grid
            %                      etc, etc. 
            gs = 2;

            k = 0.3; % length of arrows, as percentage of area of square
            nsv = beta*k*sv;

            l = 0.2;
            nad = l*k*dir2; % Note: dir2 is already scaled up by beta

            % Length of arrow head, as percentage of length of arrow staffs
            r = 0.2;

            for x = 0:gs
                for y  = 0:gs
                    xper = x/gs;
                    yper = y/gs;

                    % Arrow base
                    ab = (1-yper)*((1-xper)*p1t + xper*p2t) + yper*((1-xper)*p4t + xper*p3t);
                    % Arrow tip
                    at = point(ab.getx() + nsv.getx(), ab.gety() + nsv.gety(), ab.getz() + nsv.getz());
                    h = [h GAScene.plotline({ab, at}, 'k')];
                    am = r*ab + (1-r)*at;
                    af1 = point(am.getx() + nad.getx(), am.gety() + nad.gety(), am.getz() + nad.getz());
                    af2 = point(am.getx() - nad.getx(), am.gety() - nad.gety(), am.getz() - nad.getz());
                    h = [h GAScene.plotline({af1, at, af2}, 'k')];
                end
            end

        end

        function h = drawhairydisk(BV, c, offset)
            arguments
                BV PGA;
                c = 'b';
                offset PGA = origin(PGA);
            end

            GAScene.usefigure();
            hold on

            % TODO: perhaps remove this error message
            if ~isgrade(BV, 2)
                error('Can only draw a disk for a bivector.');
            end

            total_transform = sqrt(e12(PGA)/BV);
            if GAisa(total_transform, 'scalar') && double(e12(PGA)/BV) < 0
                f = -1;
            else
                f = 1;
            end
            
            rad = sqrt(norm(BV));
            
            
            % Creating points in a circle
            t = (0:pi/12:2*pi);
            for i=1:length(t)
                pts{i} = point(f*rad*sin(t(i)), rad*cos(t(i)), 0, PGA);
                pts{i} = inverse(total_transform) * pts{i} * total_transform;
            end
    
            h = GAScene.patch(pts, c, 0.8);
            for i=2:4:length(t)-1
                 p{1} = pts{i};
                 p{2} = 2.5*(pts{i}-pts{i-1}) + pts{i-1};
                 h = [h GAScene.plotline({p{2},p{1}},'k')];
            end
        end

        function h = drawhairyball(TV, c, offset)

            GAScene.usefigure();
            hold on
            
            v = e123coeff(TV);
            r = (abs(v)*3/4/pi)^(1./3.);
        
            [X,Y,Z] = sphere(8);

            X = r*X;
            Y = r*Y;
            Z = r*Z;


            h = plot3(X + offset.getx(), Y + offset.gety(), Z + offset.getz(), c)';
            
            for i=1:size(X,1)
                for j=1:size(X,2)
                    x(j) = X(i,j);
                    y(j) = Y(i,j);
                    z(j) = Z(i,j);
                end
                h = [h plot3(x + offset.getx(), y + offset.gety(), z + offset.getz(), c)];
            end

            if v > 0
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
                h = [h plot3(xh + offset.getx(), yh + offset.gety(), zh + offset.getz(), c)'];
            end
        end

        %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%
        %        Item Management       %
        %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%

        function manage_items(in)

            persistent static_items;
            persistent num_items;

            if isempty(num_items)
                num_items = 0;
            end

            % The user is either trying to:
            %     1. Add an item
            %     2. Delete an item
            %     3. Clear the items

            % If the user inserts a GASceneItem, they are trying to add it
            % If the user passes in a boolean, they either
            %     pass in true to clear the items
            %     pass in false to list out the items
            % If the user passes in an integer, the want to delete the item indexed by that integer

            if isa(in, "GASceneItem")
                GAScene.usefigure();
                static_items{num_items + 1} = in;
                num_items = num_items + 1;
            elseif islogical(in)
                if in
                    for i = 1:size(static_items, 2)
                        delete(static_items{i}.drawing_handles);
                    end
                    static_items = [];
                    num_items = 0;
                    
                else
                    % TODO: make scene name changable. currently GA Scene.
                    %scenename = "GA Scene";
                    %fprintf("~~~ Scene '%s' Items ~~~\n", scenename)
                    fprintf("~~~ Scene Still Items ~~~\n")
                    if isempty(static_items)
                        fprintf("   NO ITEMS TO DISPLAY\n")
                    end
                    for i = 1:size(static_items, 2)
                        fprintf("%d. %s\n", i, display_name(static_items{i}))
                    end

                    fprintf("~~~~~~~~~~~~~~~~~~~~~~~~~\n")
                    %fprintf("~~~~~~~~~~~%s~~~~~~~~~~~\n", repmat('~', 1, strlength(scenename)))

                end
            end
        end

        function displayitems()
            GAScene.manage_items(false);
            GAScene.manage_dynamic_items(false, false);
        end

        function clearitems()
            GAScene.manage_items(true);
            GAScene.manage_dynamic_items(false, true);
        end

        function additem(item)
            arguments
                item GASceneItem;
            end
            GAScene.manage_items(item);
        end

        function adddynamicitem(item)
            arguments
                item GASceneDynamicItem;
            end
            GAScene.manage_dynamic_items(false, item);
        end

        function deleteitem(itemindex)
            arguments
                itemindex uint32;
            end
            GAScene.manage_items(itemindex);
        end

        function refreshdynamicitems()
            GAScene.manage_dynamic_items(true, true);
        end

        %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%
        %    Dynamic Item Management   %
        %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%
        
        function manage_dynamic_items(is_draw, in, extra)
            arguments
                is_draw;
                in;
                extra = [];
            end
            persistent dynamic_items;
            persistent num_items;
            
            if isempty(num_items)
                num_items = 0;
            end

            if is_draw 
                % Here, we try to redraw all the items in the dynamic items list
                for i = 1:size(dynamic_items, 2)
                    dynamic_items{i} = dynamic_items{i}.update();
                end
            else
                % Here, we assume the user is trying to either
                %    1. Trying to add an item (GASceneDynamicItem)
                %    2. Trying to clear the items (logical, true)
                %    3. Trying to print the list of items (logical, false)
                %    4. Trying to delete a particular item (integer corresponding to index)

                if isa(in, "GASceneDynamicItem")
                    % Scenario 1: The user is trying to add an item
                    GAScene.usefigure();
                    dynamic_items{num_items + 1} = in;
                    num_items = num_items + 1;

                elseif islogical(in)
                    if in
                        % Scenario 2: The user is trying to clear the items
                        for i = 1:size(dynamic_items, 2)
                            delete(dynamic_items{i}.drawing_handles);
                        end
                        dynamic_items = [];
                        num_items = 0;
                        
                    else
                        % Scenario 3: The user is trying to print the list of items

                        % TODO: make scene name changable. currently GA Scene.
                        %scenename = "GA Scene";
                        %fprintf("~~~ Scene '%s' Dynamic Items ~~~\n", scenename)
                        
                        fprintf("~~~ Scene Dynamic Items ~~~\n")
                        
                        if isempty(dynamic_items)
                            fprintf("   NO ITEMS TO DISPLAY\n")
                        end
                        for i = 1:size(dynamic_items, 2)
                            fprintf("%d. %s\n", i, display_name(dynamic_items{i}))
                        end

                        fprintf("~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

                        %fprintf("~~~~~~~~~~~%s~~~~~~~~~~~\n", repmat('~', 1, strlength(scenename)))
                    end
                end

                % TODO: handle scenario 4

            end
        end


        %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%
        %       Figure Management      %
        %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%

        function obj = getfigure()
            persistent static_fig;
            
            if ~isempty(static_fig) && isvalid(static_fig)
                obj = static_fig;
            else
                GAScene.clearitems();
                % TODO: make figure name changable
                static_fig = figure('Name', "GA Scene",'NumberTitle','off');
                axis equal
                view(3)

                ax = gca;
                addlistener(ax, {'XLim', 'YLim', 'ZLim'}, 'PostSet', @(obj, evd)GAScene.manage_dynamic_items(true, obj, evd));
                
                obj = static_fig;
            end
        end

        function usefigure()
            figure(GAScene.getfigure());
            axis equal;
            view(3)
        end




        %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%
        %     ????????????????????     %
        %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%


        % TODO: Implement (convert p.m to appropriate values)
        % TODO: consider depricating(?)
        % function LineMesh(pts, c)
        %     if nargin == 1
        %         c = 'b';
        %     end
        %     M = size(pts, 1);
        %     N = size(pts, 2);
        %     for i=1:M
        %         for j=1:N
        %             p = pts{i,j};
        %             x(j) = p.m(2);
        %             y(j) = p.m(3);
        %             z(j) = p.m(4);
        %         end
        %         plot3(x,y,z,c);
        %     end
        %     for j=1:N
        %         for i=1:M
        %             p = pts{i,j};
        %             x(i) = p.m(2);
        %             y(i) = p.m(3);
        %             z(i) = p.m(4);
        %         end
        %         plot3(x,y,z,c);
        %     end
        % end

        %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%
        %        Generic Tooling       %
        %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%

        function R = boundingboxclip(xrange, yrange, zrange, p)
            arguments
                xrange;
                yrange;
                zrange;
                p;
            end
            p = p{1};

            x = p.getx();
            y = p.gety();
            z = p.getz();
            R = {point(clip(x, xrange(1), xrange(2)), clip(y, yrange(1), yrange(2)), clip(z, zrange(1), zrange(2)))};
        end

        function b = isinboundingbox(xrange, yrange, zrange, p)
            bx = xrange(1) <= p.getx() && p.getx() <= xrange(2);
            by = yrange(1) <= p.gety() && p.gety() <= yrange(2);
            bz = zrange(1) <= p.getz() && p.getz() <= zrange(2);

            b = bx && by && bz;
        end
    end

    % ******************** Private Static Methods ********************

    methods (Access = private, Static)

        %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%
        %     Private Static Tools     %
        %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%

        function a = gaarea_(px, py)
            %gaarea_(px,py): compute the area of a polygon.
            
            a = 0;
            for i = 1:length(px)-1
                a = a + px(i)*py(i+1) - px(i+1)*py(i);
            end
            a = a + px(length(px))*py(1) - px(1)*py(length(py));
            a = a/2;
        end
    end
end
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
        %           Settings           %
        %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%

        function settings()
            %SETTINGS - Displays the current configuration settings for PGABLE.
            %   To retrieve a particular setting, run "GA.[setting]".
            %   For example, to retrieve the value of autoscalar, run "GA.autoscalar".
            %   To change the value of a particular setting, run "GA.[setting]([value])".
            %   For example, to set the value of autoscalar to false, run
            %   "GA.autoscalar(false)".
            %   For more information on a particular setting, run help "GA.[setting]".
            %   To surpress the console output of changing a settings, set the second
            %   parameter to true, for example "GA.epsilon_tolerance(1E-13, true)" will set
            %   the epsilon tolerance to 1E-13 without printing the change to the console.

            disp("   ~~~~~~~~~~ Settings ~~~~~~~~~~")
            disp("   dimensions:           " + GAScene.dimensions())
            
        end

        function val = dimensions(newval, surpress_output)
            %DIMENSIONS - Set/retreive the DIMENSION setting.
            %   The DIMENSION setting is the dimensionality of the scene display.
            %   Currently supported values are 2 (for 2D) and 3 (for 3D)
            %   If no argument is provided, DIMENSIONS returns the current value of the
            %   DIMENSIONS setting.
            
            arguments
                newval = [];
                surpress_output = false;
            end

            persistent currentval;
            
            % By default the autoscalar setting is set to true
            if isempty(currentval)
                currentval = 3;
            end

            if isempty(newval)
                % User is trying to retrieve the current value
                val = currentval;
            else
                % User is trying to set the value
                if isnumeric(newval) && (newval == 2 || newval == 3)
                    currentval = newval;
                    if ~surpress_output
                        disp("   dimensions set to " + currentval)
                    end
                else
                    error('dimensions must be either 2 or 3')
                end
            end 
        end

        %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%
        %        Item Management       %
        %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%

        function displayitems()
            %DISPLAYITEMS - Displays the items of the scene in lists.

            GAScene.manage_still_items_(false);
            GAScene.manage_dynamic_items_(false, false);
        end

        function clearitems()
            %CLEARITEMS - Clears the items of the scene.
            GAScene.manage_still_items_(true);
            GAScene.manage_dynamic_items_(false, true);
        end

        function addstillitem(item)
            %ADDSTILLITEM - Adds a still item to a scene.
            %
            %   See also GASceneStillItem.

            arguments
                item GASceneStillItem;
            end
            GAScene.manage_still_items_(item);
            GAScene.refreshdynamicitems();
        end

        function adddynamicitem(item)
            %ADDDYNAMICITEM - Adds a dynamic item to a scene.
            %
            %   See also GASceneDynamicItem.

            arguments
                item GASceneDynamicItem;
            end
            GAScene.manage_dynamic_items_(false, item);
        end

        function deletestillitem(itemindex)
            %DELETESTILLITEM - Deletes a still item from a scene by the given index number.
            %
            %   See also GAScene.displayitems.

            arguments
                itemindex uint32;
            end
            GAScene.manage_still_items_(itemindex);
        end

        function deletedynamicitem(itemindex)
            %DELETEDYNAMICITEM - Deletes a dynamic item from a scene by the given index number.
            %
            %   See also GAScene.displayitems.

            arguments
                itemindex uint32;
            end
            GAScene.manage_dynamic_items_(false, itemindex);
        end

        function refreshdynamicitems()
            %REFRESHDYNAMICITEMS - A command to rerender dynamic items in the scene.

            GAScene.manage_dynamic_items_(true, true);
        end

        %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%
        %       Figure Management      %
        %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%

        function obj = getfigure()
            %GETFIGURE - Returns the figure object for GASCene.

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
                addlistener(ax, {'XLim', 'YLim', 'ZLim'}, 'PostSet', @(obj, evd)GAScene.manage_dynamic_items_(true, obj, evd));
                
                obj = static_fig;
            end
        end

        function usefigure()
            %USEFIGURE - Sets the GAScene figure in use.

            figure(GAScene.getfigure());
            axis equal;
            view(3)
        end
    end

    % ******************** Private Static Methods ********************

    methods (Access = private, Static)

        %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%
        %        Item Management       %
        %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%

        function manage_still_items_(in)

            persistent static_items;
            persistent num_items;

            if isempty(num_items)
                num_items = 0;
            end

            % The user is either trying to:
            %     1. Add an item
            %     2. Delete an item
            %     3. Clear the items

            % If the user inserts a GASceneStillItem, they are trying to add it
            % If the user passes in a boolean, they either
            %     pass in true to clear the items
            %     pass in false to list out the items
            % If the user passes in an integer, the want to delete the item indexed by that integer

            if isa(in, "GASceneStillItem")
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
            elseif isa(in, "uint32")
                delete(static_items{in}.drawing_handles);
                static_items(:, in) = [];
                num_items = num_items - 1;
                % Now display the new list to the user
                fprintf("\nItem has been Deleted. The still item list is now:\n")
                GAScene.manage_still_items_(false);
            end
        end

        function manage_dynamic_items_(is_draw, in, extra)
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
                elseif isa(in, "uint32")
                    delete(dynamic_items{in}.drawing_handles);
                    dynamic_items(:, in) = [];
                    num_items = num_items - 1;
                    % Now display the new list to the user
                    fprintf("\nItem has been Deleted. The dynamic item list is now:\n")
                    GAScene.manage_dynamic_items_(false, false);
                end

            end
        end


    end
end
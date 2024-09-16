classdef GASceneItem
    
    % TODO: Create GASceneItem, which holds a GA object and drawings associated with it.
    % TODO: Perhaps make these properties private
    properties (Access = public)
        element;  % The GA element associated with this item.
        drawing_handles; % Handles for the drawings associated with this item.
    end

    methods (Access = public)
        function obj = GASceneItem(element, drawing_handles)
            obj.element = element;
            obj.drawing_handles = drawing_handles;
        end

        function r = display_name(self)
            r = char(self.element, true);
        end

        function delete(A)
            delete(A.drawings);
        end

        
    end
end
classdef GASceneDynamicItem
    % A dynamic scene item

    % PGABLE, Copyright (c) 2024, University of Waterloo
    % Copying, use and development for non-commercial purposes permitted.
    %          All rights for commercial use reserved; for more information
    %          contact Stephen Mann (smann@uwaterloo.ca)
    %
    %          This software is unsupported.

    % TODO: Perhaps make these properties private
    % TODO: make help message of this class more useful.
    properties (Access = public)
        element;  % The GA element associated with this item.
        drawing_handles; % Handles for the drawings associated with this item.
        drawing_function;
    end

    methods (Access = public)
        function obj = GASceneDynamicItem(element, drawing_handles, drawing_function)
            obj.element = element;
            obj.drawing_handles = drawing_handles;
            obj.drawing_function = drawing_function;
        end

        function r = display_name(self)
            r = char(self.element, true);
        end

        % TODO: actually use this. perhaps use self instead of A
        function delete(A)
            delete(A.drawings);
        end

        function self = update(self)
            delete(self.drawing_handles);
            self.drawing_handles = self.drawing_function();
        end

        
    end
end
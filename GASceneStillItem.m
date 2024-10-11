classdef GASceneStillItem
    % A GA scene item.

    % PGABLE, Copyright (c) 2024, University of Waterloo
    % Copying, use and development for non-commercial purposes permitted.
    %          All rights for commercial use reserved; for more information
    %          contact Stephen Mann (smann@uwaterloo.ca)
    %
    %          This software is unsupported.
    
    % TODO: Create GASceneStillItem, which holds a GA object and drawings associated with it.
    % TODO: Perhaps make these properties private
    properties (Access = public)
        element;  % The GA element associated with this item.
        drawing_handles; % Handles for the drawings associated with this item.
    end

    methods (Access = public)
        function obj = GASceneStillItem(element, drawing_handles)
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
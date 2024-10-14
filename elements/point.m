function e = point(x, y, z, model)
    %POINT - Create a point.
    %   Either supply a (x, y, z) coordinates and the model to create a point in said
    %   model (if points are a valid geometric element in the model).
    %   Alternatively, supply an OGA vector as the first argument and the model as the
    %   second argument to create a point in the specified model whose coordinates correspond
    %   to the OGA vector.

    % PGABLE, Copyright (c) 2024, University of Waterloo
    % Copying, use and development for non-commercial purposes permitted.
    %          All rights for commercial use reserved; for more information
    %          contact Stephen Mann (smann@uwaterloo.ca)
    %
    %          This software is unsupported.
    arguments
        x;
        y = [];
        z = [];
        model = GA.model;
    end

    if nargin >= 3
        if isnumeric(x) && isnumeric(y) && isnumeric(z)
            
        else
            error("If supplying (x, y, z) to point, the values must be numeric")
        end
    elseif nargin == 1
        if isa(x, "OGA")
            z = x.getz();
            y = x.gety();
            x = x.getx();
        end
    elseif nargin == 2
        if isa(y, "GA") || isa(y, "string")
            model = y;
        end
        if isa(x, "OGA")
            z = x.getz();
            y = x.gety();
            x = x.getx();
        end
    else
       error("Incorrect number of arguments provided") 
    end

    if isa(model, "GA")
        model = modelname(model);
    end
  
    switch model
        case "PGA"
            e = PGA(0, 0, 0, [-z, y, -x, 1], 0);
        case "OGA"
            error('Points do not exist in the OGA model.')
        otherwise
            error('Cannot create element due to being in an implemented GA model.')
    end
end
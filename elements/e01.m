function e = e01(model)
    % e01

    % PGABLE, Copyright (c) 2024, University of Waterloo
    % Copying, use and development for non-commercial purposes permitted.
    %          All rights for commercial use reserved; for more information
    %          contact Stephen Mann (smann@uwaterloo.ca)
    %
    %          This software is unsupported.
    arguments
        model = GA.model;
    end
    if isa(model, "GA")
        model = modelname(model);
    end
    
    switch model
        case "PGA"
            e = PGA(0, 0, [1, 0, 0, 0, 0, 0], 0, 0);
        case "OGA"
            error('Cannot create e01 element as it does not exist in the OGA model.')
        otherwise
            error('Cannot create element due to being in an implemented GA model.')
    end
end
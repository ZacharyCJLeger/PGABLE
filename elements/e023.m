function e = e023(model)
    % e023

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
        model = model.modelname();
    end
    
    switch model
        case "PGA"
            e = PGA(0, 0, 0, [0, 0, 1, 0], 0);
        case "CGA"
            error('Cannot create e023 element as it does not exist in the CGA model.')
        case "OGA"
            error('Cannot create e023 element as it does not exist in the OGA model.')
        otherwise
            error('Cannot create element due to being in an implemented GA model.')
    end
end

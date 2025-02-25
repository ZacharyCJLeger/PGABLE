function e = e01234(model)
    % e0123

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
        case "CGA"
            e = CGA(0, 0, 0, 0, 0, 1);
        case "PGA"
            error('Cannot create e01234 element as it does not exist in the PGA model.')
        case "OGA"
            error('Cannot create e01234 element as it does not exist in the OGA model.')
        otherwise
            error('Cannot create element due to being in an implemented GA model.')
    end
end

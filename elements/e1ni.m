function e = e1ni(model)
    % e1ni

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
            e = CGA(0, 0, [0, 0, 0, 0, 0, 0, 1, 0, 0, 0], 0, 0, 0);
        case "PGA"
	    error('Cannot create e1ni element as it does not exist in the PGA model. Type  GA.model(CGA)  to switch to the CGA model.')
        case "OGA"
	    error('Cannot create e1ni element as it does not exist in the OGA model. Type  GA.model(CGA)  to switch to the CGA model.')
        otherwise
            error('Cannot create element due to being in an implemented GA model.')
    end
end

function e = origin(model)
    % origin

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
            e = no(CGA);
        case "PGA"
            e = hd(e0(PGA));
        case "OGA"
            e = OGA(0);
        otherwise
            error('Cannot create element due to being in an implemented GA model.')
    end
end

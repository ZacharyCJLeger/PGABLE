function e = origin(model)
    arguments
        model = GA.model;
    end
    if isa(model, "GA")
        model = modelname(model);
    end
  
    switch model
        case "PGA"
            e = ihd(e0(PGA));
        case "OGA"
            e = OGA(0);
        otherwise
            error('Cannot create element due to being in an implemented GA model.')
    end
end
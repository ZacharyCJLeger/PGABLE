function e = e023(model)
    arguments
        model = GA.model;
    end
    if isa(model, "GA")
        model = modelname(model);
    end
    
    switch model
        case "PGA"
            e = PGA(0, 0, 0, [0, 0, 1, 0], 0);
        case "OGA"
            error('Cannot create e023 element as it does not exist in the OGA model.')
        otherwise
            error('Cannot create element due to being in an implemented GA model.')
    end
end
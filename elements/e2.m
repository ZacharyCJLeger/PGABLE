function e = e2(model)
    arguments
        model = GA.model;
    end
    if isa(model, "GA")
        model = modelname(model);
    end
    
    switch model
        case "PGA"
            e = PGA(0, [0, 0, 1, 0], 0, 0, 0);
        case "OGA"
            e = OGA(0, [0, 1, 0], 0, 0);
        otherwise
            error('Cannot create element due to being in an implemented GA model.')
    end
end
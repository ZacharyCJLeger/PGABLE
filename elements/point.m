function e = point(x, y, z, model)
    arguments
        x (1, 1) {mustBeNumeric};
        y (1, 1) {mustBeNumeric};
        z (1, 1) {mustBeNumeric};
        model = GA.model;
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
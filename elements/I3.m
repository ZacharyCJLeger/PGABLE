function e = I3(model)
    arguments
        model = GA.model;
    end
    
    % This function returns the element e1^e2^e3, which is the pseudoscalar of OGA.

    % TODO: Perhaps throw an error when using this in a non-OGA context.
    e = e123(model);
end
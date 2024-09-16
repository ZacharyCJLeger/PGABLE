function e = I4(model)
    arguments
        model = GA.model;
    end
    
    % This function returns the element e0^e1^e2^e3, which is the pseudoscalar of PGA.

    % TODO: Perhaps throw an error when using this in a non-PGA context.
    e = e0123(model);
end
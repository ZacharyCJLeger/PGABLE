function e = I5(model)
    % I5

    % PGABLE, Copyright (c) 2025, University of Waterloo
    % Copying, use and development for non-commercial purposes permitted.
    %          All rights for commercial use reserved; for more information
    %          contact Stephen Mann (smann@uwaterloo.ca)
    %
    %          This software is unsupported.
    arguments
        model = GA.model;
    end
    
    % This function returns the element no^e1^e2^e3^ni, which is the pseudoscalar of CGA.

    % TODO: Perhaps throw an error when using this in a non-CGA context.
    e = -1*CGA.I5sign()*e12345(model);
end

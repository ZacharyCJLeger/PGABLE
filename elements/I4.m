function e = I4(model)
    % I4

    % PGABLE, Copyright (c) 2024, University of Waterloo
    % Copying, use and development for non-commercial purposes permitted.
    %          All rights for commercial use reserved; for more information
    %          contact Stephen Mann (smann@uwaterloo.ca)
    %
    %          This software is unsupported.
    arguments
        model = GA.model;
    end
    
    % This function returns the element e0^e1^e2^e3, which is the pseudoscalar of PGA.

    % TODO: Perhaps throw an error when using this in a non-PGA context.
    e = e0123(model);
end
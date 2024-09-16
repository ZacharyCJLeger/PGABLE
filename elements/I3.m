function e = I3(model)
    % I3

    % PGABLE, Copyright (c) 2024, University of Waterloo
    % Copying, use and development for non-commercial purposes permitted.
    %          All rights for commercial use reserved; for more information
    %          contact Stephen Mann (smann@uwaterloo.ca)
    %
    %          This software is unsupported.
    arguments
        model = GA.model;
    end
    
    % This function returns the element e1^e2^e3, which is the pseudoscalar of OGA.

    % TODO: Perhaps throw an error when using this in a non-OGA context.
    e = e123(model);
end
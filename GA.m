classdef (Abstract) GA
    % GA  is an abstract class that encompasses all geometric algebra objects.
    %
    %   To see settings, run GA.settings
    %
    % See also PGA, OGA.

    % PGABLE, Copyright (c) 2024, University of Waterloo
    % Copying, use and development for non-commercial purposes permitted.
    %          All rights for commercial use reserved; for more information
    %          contact Stephen Mann (smann@uwaterloo.ca)
    %
    %          This software is unsupported.

    % ******************** Public Static Methods ********************

    methods (Static = true)

        %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%
        %           Settings           %
        %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%

        function settings()
            %SETTINGS  Displays the current configuration settings for PGABLE.
            %   To retrieve a particular settings, run GA.[setting].
            %   For example, to retrieve the value of autoscalar, run GA.autoscalar.
            %   To change the value of a particular setting, run GA.[setting]([value]).
            %   For example, to set the value of autoscalar to false, run
            %   GA.autoscalar(false).
            %   For more information on a particular setting, run help GA.[setting].
            %   To surpress the console output of changing a settings, set the second
            %   parameter to true, for example GA.epsilon_tolerance(1E-13, true) will set
            %   the epsilon tolerance to 1E-13 without printing the change to the console.

            disp("   ~~~~~~~~~~ Settings ~~~~~~~~~~")
            disp("   autoscalar:           " + GA.autoscalar())
            disp("   compact_notation:     " + GA.compact_notation())
            disp("   compact_pseudoscalar: " + GA.compact_pseudoscalar())
            disp("   epsilon_tolerance:    " + GA.epsilon_tolerance())
            disp("   indicate_model:       " + GA.indicate_model())
            disp("   increasing_order:     " + GA.increasing_order())
            disp("   model:                " + GA.model())
            
        end

        % TODO: Mention in documentation that this function will now no longer convert GA elements to scalars, only scalars to GA elements
        function val = autoscalar(newval, surpress_output)
            %AUTOSCALAR  Set/retreive the AUTOSCALAR setting.
            %   The AUTOSCALAR setting is either true or false.
            %   When set to true, doubles in equations will automatically be converted to GA
            %   scalar elements.
            %   When set to false, doubles will return an error.
            %   If no argument is provided, AUTOSCALAR returns the current value of the
            %   AUTOSCALAR setting.
            
            arguments
                newval = [];
                surpress_output = false;
            end

            persistent currentval;
            
            % By default the autoscalar setting is set to true
            if isempty(currentval)
                currentval = true;
            end

            if isempty(newval)
                % User is trying to retrieve the current value
                val = currentval;
            else
                % User is trying to set the value
                if islogical(newval)
                    currentval = newval;
                    if ~surpress_output
                        disp("   autoscalar set to " + currentval)
                    end
                else
                    error('autoscalar must have value true or false.')
                end
            end 
        end

        function val = model(newval, surpress_output)
            %MODEL  Set/retreive the MODEL setting.
            %   The MODEL setting is a string indicating the current model.
            %   The argument must be either an element in the desired model or a string
            %   of the exact name of the desired model. Thus, for example, both GA.model(PGA)
            %   and GA.model("PGA") will set the current model to PGA.
            %   If no argument is provided, MODEL returns a string of the name of the
            %   current model.
            %
            %   The value of MODEL indicates which model of geometric algebra is in use.
            %   This determines, for example:
            %      - What model to construct elements such as e1
            %      - How drawing an element such as e12 should be interpreted
            %      - Whether the element e0 is a valid element or not
            %   
            %   To construct an element outside of the current model, a model can be indicated
            %   as an argument. For example, e1(PGA) will always construct e1 as a PGA
            %   element, regardless of the current selected model.
            %
            %   See also e1, origin, point.

            arguments
                newval = [];
                surpress_output = false;
            end

            persistent currentval;
            
            % By default the model is PGA
            if isempty(currentval)
                currentval = "PGA";
            end

            if isempty(newval)
                % User is trying to retrieve the current value
                val = currentval;
            else
                % User is trying to set the value

                if isa(newval, 'GA')
                    newval = modelname(newval);
                end

                currentval = newval;

                if ~surpress_output
                    disp("   model set to " + currentval)
                end
                
            end 
        end

        function val = indicate_model(newval, surpress_output)
            %INDICATE_MODEL  Set/retreive the INDICATE_MODEL setting.
            %   The INDICATE_MODEL setting is either true or false.
            %   When set to true, the current model will be displayed with the value of each
            %   GA element.
            %   When set to false, the current model will be hidden when the value of a GA
            %   element is displayed.

            arguments
                newval = [];
                surpress_output = false;
            end

            persistent currentval;
            
            % By default the indicate_model setting is set to false
            if isempty(currentval)
                currentval = false;
            end

            if isempty(newval)
                % User is trying to retrieve the current value
                val = currentval;
            else
                % User is trying to set the value
                if islogical(newval)
                    currentval = newval;
                    if ~surpress_output
                        disp("   indicate_model set to " + currentval)
                    end
                else
                    error('indicate_model must have value true or false.')
                end
            end 
        end

        function val = epsilon_tolerance(newval, surpress_output)
            %EPSILON_TOLERANCE  Set/retreive the tolerance for epsilon.
            %   The EPSILON_TOLERANCE is a non-negative real number which indicates the value
            %   for which all values whose magnitude is smaller than it will be considered
            %   epsilon. Thus, for any value x, x will be considered an epsilon if
            %                           abs(x) < GA.epsilon_tolerance.
            %   These values will be displayed as ε or -ε.
            %   By default, the value of epsilon_tolerance is 1e-15.

            arguments
                newval = [];
                surpress_output = false;
            end

            persistent currentval;
            
            % By default the indicate_model setting is set to false
            if isempty(currentval)
                currentval = 1e-15;
            end

            if isempty(newval)
                % User is trying to retrieve the current value
                val = currentval;
            else
                % User is trying to set the value
                if isnumeric(newval)
                    if newval >= 0
                        currentval = newval;
                        if ~surpress_output
                            disp("   epsilon_tolerance set to " + currentval)
                        end
                    else
                        error('epsilon_tolerance must be a non-negative number.')
                    end
                else
                    error('epsilon_tolerance must be a number.')
                end
            end 
        end

        function val = compact_notation(newval, surpress_output)
            %COMPACT_NOTATION  Set/retrieve the COMPACT_NOTATION setting.
            %   The COMPACT_NOTATION setting is either true or false.
            %   When set to true, GA elements will be displayed in compact notation.
            %   For example, the element e1^e2^e3 will be written as e123.
            %   When set to false, GA elements will be written in the full outer product form.
            %   If no argument is provided, COMPACT_NOTATION returns the current value of the
            %   COMPACT_NOTATION setting.
            %
            %   See also COMPACT_PSEUDOSCALAR.

            arguments
                newval = [];
                surpress_output = false;
            end

            persistent currentval;
            
            % By default the compact_notation setting is set to false
            if isempty(currentval)
                currentval = false;
            end

            if isempty(newval)
                % User is trying to retrieve the current value
                val = currentval;
            else
                % User is trying to set the value
                if islogical(newval)
                    currentval = newval;
                    if ~surpress_output
                        disp("   compact_notation set to " + currentval)
                    end
                else
                    error('compact_notation must have value true or false.')
                end
            end
        end

        % TODO: add comment about CGA if/when implemented
        function val = compact_pseudoscalar(newval, surpress_output)
            %COMPACT_PSEUDOSCALAR  Set/retrieve the COMPACT_PSEUDOSCALAR setting.
            %   The COMPACT_PSEUDOSCALAR setting is either true or false.
            %   When set to true, The pseudoscalar of the GA model will be represented via I[dim]
            %   where [dim] is the dimensionality of the space. Thus the following notation is used:
            %      OGA: I3 := e1^e2^e3
            %      PGA: I4 := e0^e1^e2^e3
            %   When set to false, this notation is not used.
            %   If no argument is provided, COMPACT_PSEUDOSCALAR returns the current value of the COMPACT_PSEUDOSCLAR setting.
            %
            %   See also COMPACT_NOTATION.

            arguments
                newval = [];
                surpress_output = false;
            end

            persistent currentval;
            
            % By default the compact_pseudoscalar setting is set to false
            if isempty(currentval)
                currentval = false;
            end

            if isempty(newval)
                % User is trying to retrieve the current value
                val = currentval;
            else
                % User is trying to set the value
                if islogical(newval)
                    currentval = newval;
                    if ~surpress_output
                        disp("   compact_pseudoscalar set to " + currentval)
                    end
                else
                    error('compact_pseudoscalar must have value true or false.')
                end
            end
        end
    end


    % ******************** Protected Static Methods ********************

    methods (Access = protected, Static)
        function [s_new, pl_new] = charify_val_(val, str, s, pl)
            s_new = s;
            pl_new = pl;
            
            if val ~= 0
                if val == 1
                    s_new = [s pl str];
                elseif val == -1
                    s_new = [s pl '-' str];
                elseif abs(val) < GA.epsilon_tolerance
                    if val < 0
                        s_new = [s pl '-ε*' str];
                    else
                        s_new = [s pl 'ε*' str];
                    end
                else
                    number_string = num2str(val);
                    % This converts scientific notation to use capital E to avoid confusion
                    number_string = regexprep(number_string, 'e\+?(-?\d+)', 'E$1');
                    s_new = [s pl number_string '*' str];
                end
                pl_new = ' + ';
            end
        end

        function C = getdominating_(A, B)
            if isa(A, 'GA')
                C = A;
            elseif isa(B, 'GA')
                C = B;
            else
                % TODO: return current model, perhaps. This case might not be possible.
                error('Cannot find model to cast to.')
            end
        end
    end

    % ******************** Public Methods ********************
    methods (Access = public)

        %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%
        %         Public Tools         %
        %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%

        function s = char(p, indicate_model)
            %CHAR  Returns the string representation of a GA element.
            
            arguments
                p;
                indicate_model = false;
            end
            
            if size(p, 2) > 1
                % TODO: The following below is a demo of an idea. Essentially, we could allow the user
                %       to work with matrices of GA objects rather than one at a time.
                s = [];
                for p_element = p
                    s = [s char(9) char_(p_element)];
                end
            else 
                s = char_(p);
            end

            if GA.indicate_model() || indicate_model
                s = ['(' convertStringsToChars(modelname(p)) ') ' s];
            end
        end

        function display(p)
            disp(' ');
            disp([inputname(1),' = '])
            disp(' ');
            disp(['     ' char(p)])
            disp(' ');
        end

        %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%
        %           Wrappers           %
        %%%%%%%%%%~%%%%%%%%%%~%%%%%%%%%%

        % Addition, subtraction, negation
        
        function R = plus(A, B)
            C = GA.getdominating_(A, B);
            R = plus_(C.cast(A), C.cast(B));
        end

        function R = uplus(A)
            R = A;
        end

        function R = minus(A, B)
            C = GA.getdominating_(A, B);
            R = minus_(C.cast(A), C.cast(B));
        end

        function R = uminus(A)
            R = uminus_(A);
        end

        % Outer product

        function R = outer(A, B)
            C = GA.getdominating_(A, B);
            R = outer_(C.cast(A), C.cast(B));
        end

        function R = mpower(A, B)
            C = GA.getdominating_(A, B);
            R = outer_(C.cast(A), C.cast(B));
        end

        % Inner product

        function R = inner(A, B)
            % INNER  Compute the inner product of A and B.
            %
            % See also lcont, rcont, outer, product

            C = GA.getdominating_(A, B);
            R = inner_(C.cast(A), C.cast(B));
        end

        function R = times(A, B)
            C = GA.getdominating_(A, B);
            R = inner_(C.cast(A), C.cast(B));
        end

        % Contractions

        function R = leftcontraction(A, B)
            C = GA.getdominating_(A, B);
            R = leftcontraction_(C.cast(A), C.cast(B));
        end

        function R = lcont(A, B)
            C = GA.getdominating_(A, B);
            R = leftcontraction_(C.cast(A), C.cast(B));
        end

        function R = rightcontraction(A, B)
            C = GA.getdominating_(A, B);
            R = rightcontraction_(C.cast(A), C.cast(B));
        end

        function R = rcont(A, B)
            C = GA.getdominating_(A, B);
            R = rightcontraction_(C.cast(A), C.cast(B));
        end

        % Equalities and inequalities

        function b = eq(A, B)
            C = GA.getdominating_(A, B);
            b = eq_(C.cast(A), C.cast(B));
        end

        function b = eeq(A, B)
            C = GA.getdominating_(A, B);
            b = eeq_(C.cast(A), C.cast(B));
        end

        function b = ne(A, B)
            C = GA.getdominating_(A, B);
            b = ne_(C.cast(A), C.cast(B));
        end

        % Geometric product

        function R = product(A, B)
            C = GA.getdominating_(A, B);
            R = product_(C.cast(A), C.cast(B));
        end

        function R = prod(A, B)
            C = GA.getdominating_(A, B);
            R = product_(C.cast(A), C.cast(B));
        end

        function R = mtimes(A, B)
            C = GA.getdominating_(A, B);
            R = product_(C.cast(A), C.cast(B));
        end

        % Inverse

        function R = inverse(A)
            R = inverse_(A);
        end

        % Divide

        function R = divide(A, B)
            C = GA.getdominating_(A, B);
            R = divide_(C.cast(A), C.cast(B));
        end

        function R = mrdivide(A, B)
            C = GA.getdominating_(A, B);
            R = divide_(C.cast(A), C.cast(B));
        end

        % Logs, exponentials, roots

        function R = wexp(A)
            R = wexp_(A);
        end

        function R = gexp(A)
            R = gexp_(A);
        end

        function R = glog(A)
            R = glog_(A);
        end

        function R = sqrt(A)
            R = sqrt_(A);
        end

        % Norms and normalization

        function r = norm(A)
            %NORM  returns the norm of the multivector.

            r = norm_(A);
        end

        function r = vnorm(A)
            %VNORM  Returns the vanishing norm of the multivector.

            r = vnorm_(A);
        end

        function R = normalize(A)
            %NORMALIZE  Returns the normalized multivector.

            R = normalize_(A);
        end

        % Dual

        function R = dual(A)
            %DUAL  Computes the dual.
            %
            %   See also d.

            R = dual_(A);
        end

        function R = d(A)
            %D  Shorthand for computing the dual.
            %
            %   See also dual.

            R = dual_(A);
        end

        % Inverse dual

        function R = inversedual(A)
            %INVERSEDUAL  Computes the inverse dual.
            %
            %   See also invdual, id.
            R = inversedual_(A);
        end

        function R = invdual(A)
            %INVDUAL  Shorthand for computing the inverse dual.
            %
            %   See also inversedual, id.

            R = inversedual_(A);
        end

        function R = id(A)
            %ID  Shorthand for computing the inverse dual.
            %
            %   See also inversedual, id.

            R = inversedual_(A);
        end

        % Hodge dual

        function R = hodgedual(A)
            %HODGEDUAL  Computes the hodge dual.
            %
            %   See also hdual, hd.
            R = hodgedual_(A);
        end

        function R = hdual(A)
            %HDUAL  Shorthand for computing the hodge dual.
            %
            %   See also hodgedual, hd.
            R = hodgedual_(A);
        end

        function R = hd(A)
            %HD  Shorthand for computing the hodge dual.
            %
            %   See also hodgedual, hdual.
            R = hodgedual_(A);
        end

        % Inverse Hodge dual

        function R = inversehodgedual(A)
            %INVERSEHODGEDUAL  Computes the inverse hodge dual.
            %
            %   See also invhodgedual, invhdual, ihd.

            R = inversehodgedual_(A);
        end

        function R = invhodgedual(A)
            %INVHODGEDUAL  Shorthand for computing the inverse hodge dual.
            %
            %   See also inversehodgedual, invhdual, ihd.

            R = inversehodgedual_(A);
        end

        function R = invhdual(A)
            %INVHDUAL  Shorthand for computing the inverse hodge dual.
            %
            %   See also inversehodgedual, invhodgedual, ihd.
            R = inversehodgedual_(A);
        end

        function R = ihd(A)
            %IHD  Shorthand for computing the inverse hodge dual.
            %
            %   See also inversehodgedual, invhodgedual, invhdual.
            R = inversehodgedual_(A);
        end

        % JMap

        function R = jmap(A)
            %JMAP  Computes the jmap, also called the poincare dual.
            %
            %   See also poincaredual, pdual, pd.

            R = jmap_(A);
        end

        function R = poincaredual(A)
            %POINCAREDUAL  Computes the poincare dual, also called the jmap.
            %
            %   See also jmap, pdual, pd.

            R = jmap_(A);
        end

        function R = pdual(A)
            %PDUAL  Shorthand for computing the poincare dual.
            %
            %   See also jmap, poincaredual, pd

            R = jmap_(A);
        end

        function R = pd(A)
            %PD  Shorthand for computing the poincare dual.
            %
            %   See also jmap, poincaredual, pdual.
            R = jmap_(A);
        end

        % Reverse

        function R = reverse(A)
            %REVERSE  Computes the reverse.
            %
            %   See also rev.

            R = reverse_(A);
        end

        function R = rev(A)
            %REV  Shorthand for computing the reverse.
            %
            %   See also reverse.

            R = reverse_(A);
        end

        % Meet and join

        function R = meet(A, B)
            %MEET  Computes the meet of two multivectors.

            C = GA.getdominating_(A, B);
            R = meet_(C.cast(A), C.cast(B));
        end

        function R = join(A, B)
            %JOIN  Computes the join of two multivectors.

            C = GA.getdominating_(A, B);
            R = join_(C.cast(A), C.cast(B));
        end

        % Conjugate and involution

        function R = conjugate(A)
            %CONJUGATE  Computes the conjugate.
            %
            %   See also conj.

            R = conjugate_(A);
        end

        function R = conj(A)
            %CONJ  Shorthand for computing the conjugate
            %
            %   See also conjugate.

            R = conjugate_(A);
        end

        function R = gradeinvolution(A)
            %GRADEINVOLUTION  Computing the grade involution.
            %
            %   See also gi.

            R = gradeinvolution_(A);
        end

        function R = gi(A)
            %GI  Shorthand for computing the grade involution.
            %
            %   See also gradeinvolution.

            R = gradeinvolution_(A);
        end

        % Grades

        function R = grade(A, n)
            arguments
                A GA;
                n (1, 1) int32 = -1;
            end

            R = grade_(A, n);
        end

        function b = isgrade(A, g)
            % ISGRADE  Returns true if the multivector A is a blade of grade g, and false otherwise.
            %
            %   See also PGA.GAisa, OGA.GAisa.
            arguments
                A GA;
                g (1, 1) uint32;
            end
            b = isgrade_(A, g);
        end

        % Coordinates

        function r = getx(A)
            r = getx_(A);
        end

        function r = gety(A)
            r = gety_(A);
        end

        function r = getz(A)
            r = getz_(A);
        end

        % Cleaning

        function R = zeroepsilons(A)
            %ZEROEPSILONS  Takes a multivector A and sets to zeros any 
            %   blade with coefficient size less than GA.epsilon_tolerance.
            %
            %See also GAZ.

            R = zeroepsilons_(A);
        end

        function R = GAZ(A)
            %GAZ  is shorthand for zeroepsilons.
            %
            %See also zeroepsilons.

            R = zeroepsilons_(A);
        end

        % Conversion

        function r = double(A)
            r = double_(A);
        end

        % Versor batch

        function rg = versorbatch(V, multivectors)
            arguments
                V GA;
                multivectors;
            end

            invV = inverse(V);
            rg = cell(size(multivectors, 1), size(multivectors, 2));

            for i = 1:size(multivectors, 1)
                for j = 1:size(multivectors, 2)
                    rg{i, j} = V*multivectors{i, j}*invV;
                end
            end
        end

        % TODO: proper arguments block, rename 2nd parameter, cleanup

        function rg = versorbatchiterate(V, multivector_list, iterations, include_zero)
            arguments
                V GA;
                multivector_list;
                iterations uint32 = 1;
                include_zero = false;
            end

            invV = inverse(V);
            mln = size(multivector_list, 2);

            if include_zero
                rg = cell(iterations + 1, mln);

                for mli = 1:mln
                    rg{1, mli} = multivector_list{mli};
                    for iter = 2:(iterations+1)
                        rg{iter, mli} = V * rg{iter-1, mli} * invV;
                    end
                end
            else
                rg = cell(iterations, mln);

                for mli = 1:mln
                    rg{1, mli} = V * multivector_list{mli} * invV;
                    for iter = 2:iterations
                        rg{iter, mli} = V * rg{iter-1, mli} * invV;
                    end
                end
            end
        end

    end

    % ******************** Abstract Public Methods ********************
    
    methods (Abstract, Access = public)
        %MATRIX  Returns the internal matrix representation of a multivector.
        %   This is for debugging purposes.
        matrix(A);

        %GAISA  Determines in a multivector and a string representing a type of multivector
        %   and returns true if the multivector is of that type.
        %   Which types are valid depends on the GA model. Thus, to see which types are
        %   permissable, run help [model].GAisa to see the list of options.
        %   For example, to see the options for PGA, run "help PGA.GAisa".
        %
        %   See also PGA.GAisa, OGA.GAisa.
        GAisa(A, t);

        %PGACAST  Casts a GA element directly into PGA, without geometric considerations,
        %   and removing any incompatible elements.
        PGAcast(A);

        %PGACAST  Casts a GA element directly into OGA, without geometric considerations,
        %   and removing any incompatible elements.
        OGAcast(A);

        %DRAW  Draws the GA element to the GA Scene figure.
        %   Note that the geometric interpretation of elements depends on the model.
        %   Thus, draw(e1(OGA)) may draw something different than draw(e1(PGA)).
        %   Not all elements can be drawn. You will receive an error if it cannot be drawn.
        draw(A, c);
        
        %MODELNAME  returns the name of the model of the element.
        modelname();

        %CAST  converts the input into an element of the model if an implicit conversion
        %   is possible.
        cast(A);
    end

    % ******************** Abstract Protected Methods ********************

    methods (Abstract, Access = protected)
        % TODO: Organize/order these methods

        plus_(A, B);
        minus_(A, B);
        uminus_(A);
        outer_(A, B);
        inner_(A, B);
        leftcontraction_(A, B);
        rightcontraction_(A, B);
        eq_(A, B);
        eeq_(A, B);
        ne_(A, B);
        product_(A, B);
        inverse_(A);
        divide_(A, B);
        gexp_(A);
        glog_(A);
        norm_(A);
        vnorm_(A);
        normalize_(A);
        dual_(A);
        inversedual_(A);
        hodgedual_(A);
        inversehodgedual_(A);
        jmap_(A);
        reverse_(A);
        meet_(A, B);
        join_(A, B);
        conjugate_(A);
        gradeinvolution_(A);
        grade_(A, n);
        isgrade_(A, g);
        getx_(A);
        gety_(A);
        getz_(A);
        zeroepsilons_(A);
        char_(A);
        double_(A);
        sqrt_(A);


        % TODO: Perhaps implement the methods below
        % wexp_(A);
        % GASignature
    end

end


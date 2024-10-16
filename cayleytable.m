function cayleytable(els, fn)
    %CAYLEYTABLE - When given an array of elements and a binary function, displays a Cayley
    %   table.

    tabl = [];
    for i = 1:length(els)
        col = [];
        for j = 1:length(els)
            res = fn(els(j), els(i));
            col = [col convertCharsToStrings(char(res))];
        end
        tabl = [tabl; pad(col)];
    end

    for i = 1:length(els)
        for j = 1:length(els)
            fprintf('%s  ', tabl(j, i));
        end
        fprintf('\n')
    end
    
end